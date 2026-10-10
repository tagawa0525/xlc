//! リポジトリの文書の整合性の検査（#0013）。
//!
//! git が追跡しているファイルを一度だけ読み、xlc_docs::Docs の検査に違反がないことを確かめる。閉じた issue の
//! `resolution: コミット <hash>` は、git で HEAD（マージの途中ではマージするブランチも）の祖先のコミットを指すことを
//! 確かめる。違反はすべて並べて示す。

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

impl Repo<'_> {
    /// リポジトリの根で git を実行する。git がなければエラーにする。
    fn git(&self, args: &[&str]) -> Result<Output> {
        let mut command = Command::new("git");
        command.arg("-C").arg(self.root).args(args);
        if self.isolated {
            for (key, _) in std::env::vars_os() {
                if key.to_string_lossy().starts_with("GIT_") {
                    command.env_remove(key);
                }
            }
        }
        command
            .output()
            .map_err(|e| format!("cannot run git {args:?} (is git installed?): {e}").into())
    }
}

/// git が追跡しているファイル（`/` 区切りのパス）。
fn tracked_files(repo: &Repo) -> Result<Vec<String>> {
    let output = repo.git(&["ls-files", "-z"])?;
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

/// `resolution` のコミットがその祖先であるべき先頭。ふつうは HEAD だけ。マージの途中（MERGE_HEAD がある）では、
/// マージするブランチも加える: マージのコミットを作る前の検査（hook）では HEAD がまだ元のブランチで、issue を
/// 閉じたブランチの上のコミットは MERGE_HEAD の祖先にしかない。衝突で止まったマージの途中でも同じく認めるが、
/// そのブランチは HEAD の祖先になる前提で、マージをやめて rebase したら、次の検査で古い hash が見つかる。
/// MERGE_HEAD が複数行になるタコ足のマージでは、最初のブランチだけを見る（このリポジトリでは使わない）
fn heads(repo: &Repo) -> Result<Vec<&'static str>> {
    let merging = repo
        .git(&["rev-parse", "--quiet", "--verify", "MERGE_HEAD"])?
        .status
        .success();
    Ok(if merging {
        vec!["HEAD", "MERGE_HEAD"]
    } else {
        vec!["HEAD"]
    })
}

/// `hash` が `heads` のどれかの祖先のコミットを指すか。違えば理由を返す。
fn commit_problem(repo: &Repo, hash: &str, heads: &[&str]) -> Result<Option<String>> {
    let commit = format!("{hash}^{{commit}}");
    if !repo
        .git(&["rev-parse", "--quiet", "--verify", &commit])?
        .status
        .success()
    {
        return Ok(Some(format!("commit {hash} does not exist")));
    }
    for head in heads {
        if is_ancestor(repo, hash, head)? {
            return Ok(None);
        }
    }
    Ok(Some(format!(
        "commit {hash} is not an ancestor of {} (rebased? update the hash)",
        heads.join(" or ")
    )))
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

#[test]
fn docs_are_consistent() -> Result<()> {
    let root = repo_root();
    let repo = Repo {
        root: &root,
        isolated: false,
    };
    let mut texts = BTreeMap::new();
    let mut others = Vec::new();
    for path in tracked_files(&repo)? {
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
    let heads = heads(&repo)?;
    for (path, hash) in docs.resolution_commits() {
        if let Some(problem) = commit_problem(&repo, &hash, &heads)? {
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
        let repo = self.repo();
        commit_problem(&repo, hash, &heads(&repo)?)
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
