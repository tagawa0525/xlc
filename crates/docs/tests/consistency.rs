//! リポジトリの文書の整合性の検査（#0013）。
//!
//! git が追跡しているファイルを一度だけ読み、xlc_docs::Docs の検査に違反がないことを確かめる。閉じた issue の
//! `resolution: コミット <hash>` は、git で HEAD かローカルのどれかのブランチの祖先のコミットを指すことを確かめる。
//! 違反はすべて並べて示す。

// テストでは assert! 系を使う（本番のコードでは clippy.toml で禁止している）
#![cfg_attr(test, allow(clippy::disallowed_macros))]

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use xlc_docs::Docs;

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
    /// リポジトリの根で git を実行する。git がなければエラーにする。
    fn git(&self, args: &[&str]) -> Result<Output> {
        let mut command = Command::new("git");
        command.arg("-C").arg(self.root).args(args);
        if self.isolated {
            for key in inherited_git_variables() {
                command.env_remove(key);
            }
        }
        command
            .output()
            .map_err(|e| format!("cannot run git {args:?} (is git installed?): {e}").into())
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

/// `hash` が HEAD か、ローカルのどれかのブランチの祖先のコミットを指すか。違えば理由を返す。マージのコミットを作る
/// 前の検査（pre-merge-commit の hook）では、HEAD はまだ元のブランチで、MERGE_HEAD もまだない（git 2.55 で確かめた）。
/// issue を閉じたブランチの上のコミットは、そのブランチの祖先として認める。rebase で捨てた古いコミットは、どの
/// ブランチの祖先でもなくなるので見つかる。マージしていないブランチの上のコミットは、そのブランチを消すまで見逃す
fn commit_problem(repo: &Repo, hash: &str) -> Result<Option<String>> {
    let commit = format!("{hash}^{{commit}}");
    if !repo
        .git(&["rev-parse", "--quiet", "--verify", &commit])?
        .status
        .success()
    {
        return Ok(Some(format!("commit {hash} does not exist")));
    }
    if is_ancestor(repo, hash, "HEAD")? {
        return Ok(None);
    }
    let branches = repo.git(&[
        "for-each-ref",
        "--contains",
        hash,
        "--format=%(refname)",
        "refs/heads/",
    ])?;
    if !branches.status.success() {
        return Err(format!(
            "git for-each-ref failed: {}",
            String::from_utf8_lossy(&branches.stderr)
        )
        .into());
    }
    Ok(if branches.stdout.is_empty() {
        Some(format!(
            "commit {hash} is not an ancestor of HEAD or any branch (rebased? update the hash)"
        ))
    } else {
        None
    })
}

/// `hash` が `head` の祖先か。
fn is_ancestor(repo: &Repo, hash: &str, head: &str) -> Result<bool> {
    let ancestor = repo.git(&["merge-base", "--is-ancestor", hash, head])?;
    match ancestor.status.code() {
        Some(0) => Ok(true),
        Some(1) => Ok(false),
        _ => Err(format!(
            "git merge-base failed: {}",
            String::from_utf8_lossy(&ancestor.stderr)
        )
        .into()),
    }
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
    for (path, hash) in docs.resolution_commits() {
        if let Some(problem) = commit_problem(&repo, &hash)? {
            violations.push(format!("{path}: `resolution` names {problem}"));
        }
    }
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

/// 検証の記録の表が、比較の段階の JSON から描いた表と同じ（ADR-0176）。JSON は DVC の出力なので、無い worktree では
/// `dvc checkout examples/<問題>/dvc.yaml` で戻す（落ちたときのメッセージに出る）。
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

    fn problem(&self, hash: &str) -> Result<Option<String>> {
        commit_problem(&self.repo(), hash)
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
/// ないと DVC が worktree の中に別のキャッシュを黙って作り、ENDF を deps にする段階が止まる（ADR-0123、ADR-0126）。
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
