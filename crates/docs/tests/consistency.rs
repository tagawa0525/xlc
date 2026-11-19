//! リポジトリの文書の整合性の検査（#0013）。
//!
//! git が追跡しているファイルを一度だけ読み、xlc_docs::Docs の検査に違反がないことを確かめる。閉じた issue の
//! `resolution: コミット <hash>` と、文書が引くコミットの hash（`Docs::cited_hashes` のうちコミットを指すもの）は、
//! git で HEAD、ローカルのブランチ、タグのどれかから辿れることを確かめる。違反はすべて並べて示す。

// テストでは assert! 系を使う（本番のコードでは clippy.toml で禁止している）
#![cfg_attr(test, allow(clippy::disallowed_macros))]

use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

use xlc_docs::{Citation, Docs};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

/// 中身を読むファイルの拡張子（文書と、参照を書くソースとスクリプト）。
const TEXT_EXTENSIONS: [&str; 5] = ["md", "rs", "py", "toml", "sh"];

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..").join("..")
}

/// git を実行するリポジトリ。
struct Repo<'a> {
    root: &'a Path,
    /// 親の git（hook）から受け継いだ GIT_ で始まる環境変数（GIT_DIR、GIT_INDEX_FILE など）を外す。テストで作る
    /// 一時的なリポジトリに使う。外さないと、hook の中で走るテストの操作が、hook を呼んだ本物のリポジトリに効く。
    /// 本物のリポジトリでは外さない（hook の中では、GIT_INDEX_FILE がマージの結果の index を指す）
    isolated: bool,
}

/// 呼び出し元の環境にある `GIT_*` の名前。
#[allow(
    clippy::disallowed_methods,
    reason = "テスト用のリポジトリに、呼び出し元の GIT_* を持ち込まない（ADR-0137）"
)]
fn inherited_git_variables() -> Vec<std::ffi::OsString> {
    std::env::vars_os()
        .map(|(key, _)| key)
        .filter(|key| key.to_string_lossy().starts_with("GIT_"))
        .collect()
}

impl Repo<'_> {
    /// リポジトリの根で git を実行するコマンド。
    fn command(&self, args: &[&str]) -> Command {
        let mut command = Command::new("git");
        command.arg("-C").arg(self.root).args(args);
        if self.isolated {
            for key in inherited_git_variables() {
                command.env_remove(key);
            }
        }
        command
    }

    /// リポジトリの根で git を実行する。git がなければエラーにする。
    fn git(&self, args: &[&str]) -> Result<Output> {
        self.command(args)
            .output()
            .map_err(|e| format!("cannot run git {args:?} (is git installed?): {e}").into())
    }

    /// 標準入力に `input` を渡して git を実行し、成功したときの標準出力を返す。入力は別のスレッドで書く（git が
    /// 書く出力でパイプが詰まると、書き終わるまで待つ側と読まれるのを待つ側で止まる）
    fn stdout_with_input(&self, args: &[&str], input: String) -> Result<String> {
        let mut child = self
            .command(args)
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .map_err(|e| format!("cannot run git {args:?} (is git installed?): {e}"))?;
        let mut stdin = child.stdin.take().ok_or("git has no stdin")?;
        let writer = std::thread::spawn(move || stdin.write_all(input.as_bytes()));
        let output = child.wait_with_output()?;
        let written = writer
            .join()
            .map_err(|_| "the thread writing to git panicked")?;
        // git が先に止まると書き込みは broken pipe になるので、git の失敗（標準エラーに原因）を先に示す
        if !output.status.success() {
            return Err(format!(
                "git {args:?} failed: {}",
                String::from_utf8_lossy(&output.stderr)
            )
            .into());
        }
        written?;
        Ok(String::from_utf8(output.stdout)?)
    }

    /// git を実行し、成功したときの標準出力を返す。
    fn stdout(&self, args: &[&str]) -> Result<String> {
        let output = self.git(args)?;
        if !output.status.success() {
            return Err(format!(
                "git {args:?} failed: {}",
                String::from_utf8_lossy(&output.stderr)
            )
            .into());
        }
        Ok(String::from_utf8(output.stdout)?)
    }
}

/// git が追跡しているファイルのうち `pathspecs` に当たるもの（`/` 区切りのパス。空ならすべて）。
fn tracked_files(repo: &Repo, pathspecs: &[&str]) -> Result<Vec<String>> {
    let mut args = vec!["ls-files", "-z", "--"];
    args.extend_from_slice(pathspecs);
    let output = repo.git(&args)?;
    if !output.status.success() {
        return Err(format!(
            "git ls-files failed: {}",
            String::from_utf8_lossy(&output.stderr)
        )
        .into());
    }
    Ok(String::from_utf8(output.stdout)?
        .split('\0')
        .filter(|p| !p.is_empty())
        .map(str::to_string)
        .collect())
}

/// HEAD、ローカルのブランチ、タグから辿れるコミット（完全な hash）。マージのコミットを作る前の検査
/// （pre-merge-commit の hook）では、HEAD はまだ元のブランチで、MERGE_HEAD もまだない（git 2.55 で確かめた）ので、
/// マージするブランチの上のコミットは、そのブランチから辿れるものとして認める。rebase で捨てた古いコミットは、
/// どの ref からも辿れなくなるので見つかる。マージしていないブランチの上のコミットは、そのブランチを消すまで見逃す。
/// タグは、ブランチを消した後も残すと決めたコミット（`record/`、`archive/`）を辿れるようにする。reflog は数えない
/// （git gc が期限で消し、そのあとは hash を checkout できない。ADR-0122）
fn reachable_commits(repo: &Repo) -> Result<HashSet<String>> {
    Ok(repo
        .stdout(&["rev-list", "HEAD", "--branches", "--tags"])?
        .lines()
        .map(str::to_string)
        .collect())
}

/// 文書の中の hash が指すもの。
enum Named {
    /// コミット（完全な hash）。注釈つきのタグのオブジェクトは、その先のコミット
    Commit(String),
    /// 先頭が同じコミットが 2 つ以上ある
    Ambiguous,
    /// コミットでない（md5、ほかのリポジトリのコミット、blob、消えたコミット）
    Other,
}

/// `hashes` が指すものを、1 つの `git cat-file --batch-check` で引く。`<hash>^{commit}` で引くので、git は
/// コミットだけで短い hash を決める（blob や tree と先頭が重なっても ambiguous にならない）
fn name_commits<'a>(repo: &Repo, hashes: &BTreeSet<&'a str>) -> Result<BTreeMap<&'a str, Named>> {
    if hashes.is_empty() {
        return Ok(BTreeMap::new());
    }
    let input: String = hashes.iter().map(|h| format!("{h}^{{commit}}\n")).collect();
    let output = repo.stdout_with_input(
        &[
            "cat-file",
            "--batch-check=%(objectname) %(objecttype)",
            "--buffer",
        ],
        input,
    )?;
    let lines: Vec<&str> = output.lines().collect();
    if lines.len() != hashes.len() {
        return Err(format!(
            "git cat-file answered {} lines for {} names",
            lines.len(),
            hashes.len()
        )
        .into());
    }
    Ok(hashes
        .iter()
        .zip(lines)
        .map(|(hash, line)| {
            let named = match line.split_once(' ') {
                Some((full, "commit")) => Named::Commit(full.to_string()),
                Some((_, "ambiguous")) => Named::Ambiguous,
                _ => Named::Other,
            };
            (*hash, named)
        })
        .collect())
}

/// 閉じた issue の `resolution` の hash（パスと hash）がコミットを指すか、文書が引く hash（`cited`）のうち
/// コミットを指すものが HEAD、ブランチ、タグから辿れるかを確かめ、違反の文を返す。`resolution` の行も
/// `cited` に入るので、辿れるかは `cited` の側だけで見る（同じ hash を 2 度示さない）。コミットを指さない
/// 16 進は md5 などと区別できないので、`resolution` のほかでは見ない（git gc の後の消えたコミットも見逃す）
fn commit_problems(
    repo: &Repo,
    resolutions: &[(String, String)],
    cited: &[Citation],
) -> Result<Vec<String>> {
    let hashes: BTreeSet<&str> = resolutions
        .iter()
        .map(|(_, h)| h.as_str())
        .chain(cited.iter().map(|c| c.hash.as_str()))
        .collect();
    let named = name_commits(repo, &hashes)?;
    let reachable = reachable_commits(repo)?;
    let mut problems: Vec<String> = resolutions
        .iter()
        .filter(|(_, hash)| !matches!(named.get(hash.as_str()), Some(Named::Commit(_))))
        .map(|(path, hash)| format!("{path}: `resolution` names commit {hash} that does not exist"))
        .collect();
    for c in cited {
        let problem = match named.get(c.hash.as_str()) {
            Some(Named::Commit(full)) if !reachable.contains(full) => format!(
                "commit {} is not reachable from HEAD, a branch or a tag \
                 (cite the rebased commit, or tag the commit to keep it)",
                c.hash
            ),
            Some(Named::Ambiguous) => {
                format!("{} names more than one commit (write more digits)", c.hash)
            }
            _ => continue,
        };
        problems.push(format!("{}:{}: {problem}", c.path, c.line));
    }
    Ok(problems)
}

/// 今のブランチが main のとき、`trials/` の下で追跡しているファイル。ほかのブランチと detached HEAD では空。
/// ブランチの上では、記録に数値を書いた試行の入力を追跡してよく、マージの前に消す（ADR-0122）。自動のマージで
/// コミットを作る前の hook では、HEAD はまだ main で、.git/index がマージの結果なので（git 2.55 で確かめた）、
/// `git ls-files` で残りが見つかる。競合を直して `git commit -a` で締めるときに git が hook に渡す一時的な index は、
/// run-checks.sh が GIT_INDEX_FILE を外すので読まない
fn trials_tracked_on_main(repo: &Repo) -> Result<Vec<String>> {
    let branch = repo.git(&["branch", "--show-current"])?;
    if !branch.status.success() {
        return Err(format!(
            "git branch --show-current failed: {}",
            String::from_utf8_lossy(&branch.stderr)
        )
        .into());
    }
    if String::from_utf8(branch.stdout)?.trim() != "main" {
        return Ok(Vec::new());
    }
    tracked_files(repo, &["trials/"])
}

#[test]
fn docs_are_consistent() -> Result<()> {
    let root = repo_root();
    let repo = Repo {
        root: &root,
        isolated: false,
    };
    let mut texts = BTreeMap::new();
    let mut others = Vec::new();
    for path in tracked_files(&repo, &[])? {
        let extension = path.rsplit_once('.').map_or("", |(_, e)| e);
        // 追跡していても作業ツリーで消したファイルは、読まずに数えない
        let file = root.join(&path);
        if !file.exists() {
            continue;
        }
        if TEXT_EXTENSIONS.contains(&extension) {
            texts.insert(path, fs::read_to_string(file)?);
        } else {
            others.push(path);
        }
    }
    assert!(
        texts.contains_key("CLAUDE.md") && texts.contains_key("docs/adr/README.md"),
        "the repository root was not found at {}",
        root.display()
    );

    let docs = Docs::new(texts, others);
    let mut violations: Vec<String> = docs.check().iter().map(ToString::to_string).collect();
    violations.extend(commit_problems(
        &repo,
        &docs.resolution_commits(),
        &docs.cited_hashes(),
    )?);
    assert!(
        violations.is_empty(),
        "{} documentation violation(s):\n{}",
        violations.len(),
        violations.join("\n")
    );
    Ok(())
}

#[test]
fn main_tracks_no_trials() -> Result<()> {
    let root = repo_root();
    let repo = Repo {
        root: &root,
        isolated: false,
    };
    let trials = trials_tracked_on_main(&repo)?;
    assert!(
        trials.is_empty(),
        "main tracks trials (remove them, or git mv the ones to keep into examples/, before the merge; ADR-0122):\n{}",
        trials.join("\n")
    );
    Ok(())
}

/// 本番のファイル（`crates/*/src` の `tests.rs` 以外）にテストのコードがない（ADR-0143、#0089）。
#[test]
fn production_sources_hold_no_test_code() -> Result<()> {
    let root = repo_root();
    let output = Command::new("bash")
        .arg(root.join("scripts/check-no-inline-tests.sh"))
        .arg(&root)
        .output()
        .map_err(|e| format!("cannot run scripts/check-no-inline-tests.sh: {e}"))?;
    assert!(
        output.status.success(),
        "{}{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    Ok(())
}

/// 検証の記録の表が、比較のステージの JSON から描いた表と同じ（ADR-0176）。JSON は数 KB の要約で、DVC の出力
/// （`cache: false`）のまま git で追跡するので、worktree でも nix の sandbox でも `dvc checkout` なしに走る。
#[test]
fn record_tables_match_the_compare_outputs() -> Result<()> {
    let root = repo_root();
    let output = Command::new("python3")
        .arg(root.join("examples/verification_tables.py"))
        .arg("--check")
        .output()
        .map_err(|e| {
            format!("cannot run examples/verification_tables.py (is python3 installed?): {e}")
        })?;
    assert!(
        output.status.success(),
        "{}{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    Ok(())
}

/// テストで作る一時的なリポジトリ。落ちたときも含めて、片付けで消す。
struct TempRepo {
    dir: PathBuf,
}

impl TempRepo {
    /// プロセスと時刻で名前を一意にした、空のリポジトリ
    fn new(name: &str) -> Result<Self> {
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)?
            .as_nanos();
        let dir =
            std::env::temp_dir().join(format!("xlc-docs-{name}-{}-{nanos}", std::process::id()));
        fs::create_dir(&dir)?;
        let repo = TempRepo { dir };
        repo.run(&["init", "-q", "-b", "main"])?;
        Ok(repo)
    }

    fn repo(&self) -> Repo<'_> {
        Repo {
            root: &self.dir,
            isolated: true,
        }
    }

    /// `args` を実行し、成功を確かめて標準出力を返す。署名の設定は使わず、利用者の hook は `--no-verify` で
    /// 飛ばす（呼ぶ側で付ける）
    fn run(&self, args: &[&str]) -> Result<String> {
        let mut all = vec![
            "-c",
            "user.name=t",
            "-c",
            "user.email=t@t",
            "-c",
            "commit.gpgsign=false",
        ];
        all.extend_from_slice(args);
        let output = self.repo().git(&all)?;
        if !output.status.success() {
            return Err(
                format!("git {args:?}: {}", String::from_utf8_lossy(&output.stderr)).into(),
            );
        }
        Ok(String::from_utf8(output.stdout)?.trim().to_string())
    }

    /// 空のコミットを作り、その hash を返す
    fn commit(&self, message: &str) -> Result<String> {
        self.run(&[
            "commit",
            "-q",
            "--no-verify",
            "--allow-empty",
            "-m",
            message,
        ])?;
        self.run(&["rev-parse", "HEAD"])
    }

    /// 文書が完全な `hash` を引いたときの最初の違反
    fn problem(&self, hash: &str) -> Result<Option<String>> {
        Ok(self.cited(&format!("{hash}\n"))?.into_iter().next())
    }

    /// 1 つの文書 `text` が引く hash の違反
    fn cited(&self, text: &str) -> Result<Vec<String>> {
        let docs = Docs::new(
            BTreeMap::from([("docs/a.md".to_string(), text.to_string())]),
            Vec::<String>::new(),
        );
        commit_problems(&self.repo(), &[], &docs.cited_hashes())
    }
}

impl Drop for TempRepo {
    fn drop(&mut self) {
        // 片付けに失敗しても、テストの結果は変えない（一時的なディレクトリが残るだけ）
        let _ = fs::remove_dir_all(&self.dir);
    }
}

#[test]
fn branch_commits_count_until_the_branch_is_gone() -> Result<()> {
    // issue を閉じたブランチの上のコミットは、マージの前（hook の検査）もあとも認める。rebase で捨てたコミット
    // （どのブランチにもない）は認めない
    let t = TempRepo::new("branches")?;
    t.commit("chore: base")?;
    t.run(&["switch", "-q", "-c", "topic"])?;
    let fix = t.commit("fix: topic")?;
    t.run(&["switch", "-q", "-c", "gone"])?;
    let stale = t.commit("fix: gone")?;
    t.run(&["switch", "-q", "main"])?;
    t.run(&["branch", "-q", "-D", "gone"])?;
    t.commit("chore: other")?;
    assert_eq!(t.problem(&fix)?, None, "a commit on an unmerged branch");
    assert!(t.problem(&stale)?.is_some(), "a commit on no branch");
    t.run(&[
        "merge",
        "-q",
        "--no-ff",
        "--no-verify",
        "-m",
        "chore: merge",
        "topic",
    ])?;
    t.run(&["branch", "-q", "-d", "topic"])?;
    assert_eq!(
        t.problem(&fix)?,
        None,
        "a merged commit after the branch is deleted"
    );
    assert!(
        t.problem(&stale)?.is_some(),
        "a commit on no branch after the merge"
    );
    Ok(())
}

#[test]
fn tagged_commits_count() -> Result<()> {
    // ブランチを消しても、タグ（record/ などで残すと決めたコミット）から辿れるコミットは認める
    let t = TempRepo::new("tags")?;
    t.commit("chore: base")?;
    t.run(&["switch", "-q", "-c", "kept"])?;
    let kept = t.commit("fix: kept")?;
    t.run(&["tag", "-a", "-m", "keep", "record/kept", &kept])?;
    t.run(&["switch", "-q", "main"])?;
    t.run(&["branch", "-q", "-D", "kept"])?;
    assert_eq!(t.problem(&kept)?, None);
    Ok(())
}

#[test]
fn cited_commits_must_be_reachable() -> Result<()> {
    // 文書が引くコミットは、HEAD、ブランチ、タグのどれかから辿れる。コミットでない 16 進（md5、ほかの
    // リポジトリのコミット、blob）は見ない
    let t = TempRepo::new("cited")?;
    let base = t.commit("chore: base")?;
    t.run(&["switch", "-q", "-c", "gone"])?;
    let stale = t.commit("fix: gone")?;
    t.run(&["switch", "-q", "-c", "topic", "main"])?;
    let topic = t.commit("fix: topic")?;
    t.run(&["switch", "-q", "-c", "kept"])?;
    let kept = t.commit("fix: kept")?;
    t.run(&["tag", "-a", "-m", "keep", "record/kept", &kept])?;
    t.run(&["switch", "-q", "main"])?;
    t.run(&["branch", "-q", "-D", "gone", "kept"])?;
    let blob = t.run(&["hash-object", "-w", "--stdin"])?;
    let md5 = "0123456789abcdef0123456789abcdef";
    let text = format!(
        "xlc {}（分岐点は main の {}）、{}、md5 {md5}、{}、{}\n",
        &topic[..9],
        &base[..9],
        &kept[..9],
        "abcdef12",
        &blob[..9]
    );
    assert_eq!(t.cited(&text)?, Vec::<String>::new());
    let problems = t.cited(&format!("{text}rebase の前の {}\n", &stale[..9]))?;
    assert_eq!(problems.len(), 1, "{problems:?}");
    assert!(
        problems[0].starts_with("docs/a.md:2: ") && problems[0].contains(&stale[..9]),
        "{problems:?}"
    );
    Ok(())
}

#[test]
fn resolutions_name_existing_commits_and_report_once() -> Result<()> {
    // `resolution` の hash は、コミットが無ければ違反。辿れないコミットは、文書の行として 1 度だけ示す
    let t = TempRepo::new("resolution")?;
    t.commit("chore: base")?;
    t.run(&["switch", "-q", "-c", "gone"])?;
    let stale = t.commit("fix: gone")?;
    t.run(&["switch", "-q", "main"])?;
    t.run(&["branch", "-q", "-D", "gone"])?;
    let path = "docs/issues/0001-a.md".to_string();
    let text = format!("---\nresolution: コミット {}\n---\n", &stale[..9]);
    let docs = Docs::new(BTreeMap::from([(path.clone(), text)]), Vec::<String>::new());
    let problems = commit_problems(
        &t.repo(),
        &[(path.clone(), stale[..9].to_string())],
        &docs.cited_hashes(),
    )?;
    assert_eq!(problems.len(), 1, "{problems:?}");
    assert!(
        problems[0].starts_with(&format!("{path}:2: commit ")),
        "{problems:?}"
    );
    let missing = "abcdef1234";
    assert_eq!(
        commit_problems(&t.repo(), &[(path.clone(), missing.to_string())], &[])?,
        [format!(
            "{path}: `resolution` names commit {missing} that does not exist"
        )]
    );
    Ok(())
}

#[test]
fn prefixes_of_two_commits_are_ambiguous() -> Result<()> {
    // 7 桁で重なるコミットは作りにくいので、name_commits に 4 桁の先頭を渡す。日付を固めた fast-import で
    // 作るので、コミットの hash（と重なる先頭）は毎回同じ
    let t = TempRepo::new("ambiguous")?;
    let stream: String = (0..1000)
        .map(|i| {
            format!(
                "commit refs/heads/main\ncommitter t <t@t> 1700000000 +0000\ndata {}\nc{i:04}\n\n",
                "c0000".len() + 1
            )
        })
        .collect();
    t.repo()
        .stdout_with_input(&["fast-import", "--quiet"], stream)?;
    let commits = t.run(&["rev-list", "main"])?;
    let mut seen = BTreeMap::new();
    let prefix = commits
        .lines()
        .map(|c| &c[..4])
        .find(|p| seen.insert(*p, ()).is_some())
        .ok_or("no two of the commits share 4 leading digits")?;
    let one = commits.lines().next().ok_or("no commit")?;
    let named = name_commits(&t.repo(), &BTreeSet::from([prefix, &one[..9]]))?;
    assert!(matches!(named[prefix], Named::Ambiguous));
    assert!(matches!(&named[&one[..9]], Named::Commit(full) if full == one));
    Ok(())
}

#[test]
fn main_reports_trials_left_by_a_merge() -> Result<()> {
    // ブランチの上では、記録に書いた試行を trials/ に追跡してよい。マージのコミットを作る前（index がマージの
    // 結果で、HEAD はまだ main）に trials/ が残っていれば見つける
    let t = TempRepo::new("trials")?;
    t.commit("chore: base")?;
    assert_eq!(trials_tracked_on_main(&t.repo())?, Vec::<String>::new());
    t.run(&["switch", "-q", "-c", "topic"])?;
    fs::create_dir(t.dir.join("trials"))?;
    fs::write(t.dir.join("trials/input.toml"), "")?;
    t.run(&["add", "trials/input.toml"])?;
    t.commit("docs: record a trial")?;
    assert_eq!(
        trials_tracked_on_main(&t.repo())?,
        Vec::<String>::new(),
        "a trial on a topic branch"
    );
    t.run(&["switch", "-q", "main"])?;
    t.run(&["merge", "-q", "--no-ff", "--no-commit", "topic"])?;
    assert_eq!(
        trials_tracked_on_main(&t.repo())?,
        vec!["trials/input.toml".to_string()],
        "a trial in the merge result on main"
    );
    Ok(())
}

/// linked worktree の `references` は、main の worktree の `references/` を指す symlink でなければならない。
/// ないと DVC が worktree の中に別のキャッシュを黙って作り、ENDF を deps にするステージが止まる（ADR-0123、ADR-0126）。
/// main の worktree では `references/` は実体なので、確かめることはない。
#[test]
fn linked_worktrees_share_references() -> Result<()> {
    let root = repo_root().canonicalize()?;
    let repo = Repo {
        root: &root,
        isolated: false,
    };
    // main の worktree は一覧の先頭（worktree-add.sh と同じ求め方）
    let list = repo.stdout(&["worktree", "list", "--porcelain"])?;
    let main_root = list
        .lines()
        .next()
        .and_then(|l| l.strip_prefix("worktree "))
        .ok_or("git worktree list names no worktree")?;
    let main_root = PathBuf::from(main_root).canonicalize()?;
    if main_root == root {
        return Ok(());
    }
    let main_references = main_root.join("references");
    assert!(
        main_references.is_dir(),
        "{} does not exist; build it in the main worktree first (docs/environment.md; ADR-0124)",
        main_references.display()
    );
    let references = root.join("references");
    let hint = "link it from the worktree root: \
                ln -s \"$(git rev-parse --path-format=absolute --git-common-dir)/../references\" references";
    match fs::symlink_metadata(&references) {
        Err(_) => panic!("{} is missing (ADR-0126); {hint}", references.display()),
        Ok(m) => assert!(
            m.file_type().is_symlink(),
            "{} exists but is not a symlink to the main worktree's references/ (ADR-0126); \
             move it away, then {hint}",
            references.display()
        ),
    }
    assert_eq!(
        fs::canonicalize(&references)
            .map_err(|e| format!("{} does not resolve ({e}); {hint}", references.display()))?,
        main_references.canonicalize()?,
        "{} does not resolve to {} (ADR-0126); remove it, then {hint}",
        references.display(),
        main_references.display()
    );
    Ok(())
}
