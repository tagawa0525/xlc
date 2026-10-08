# xlc

PWR の解析コード。xs（ENDF から多群定数）→ lattice（MOC の格子計算）→ core（3次元の NEM）の段階からなる。

## 文書

- 決定と理由: docs/adr/（書き方と運用は docs/adr/README.md）
- 段階の目的・範囲・作業の順序・検証の対象: docs/plans/
- 未解決の問題とやること: docs/issues/（書き方と運用は docs/issues/README.md。ブランチ名とコミットの `Refs:`、
  ADR との両方向の対応を含む）
- 検証の記録: examples/*/verification.md。時系列に節を足し、過去の節は書き換えない。誤りは日付つきの注記で直す。
  各節に、測ったときの xlc のコミットと日付を書く（コミットが分からない過去の測定は日付だけ）
- 検証の環境（NJOY、OpenMC、核データ）: docs/environment.md。`references/` は git で追跡しない

## 作業の流れ

- GitHub のリモートはない。PR の代わりに /git-merge でマージする（rebase、チェック、/code-review と指摘への対応、
  `git merge --no-ff`、後片付けを含む）
- チェックは `~/.claude/skills/language-checks/scripts/run-checks.sh` で、終了コードで判定する。`checks | tail && git merge`
  のようにつなぐと、失敗してもマージに進む
- worktree は `~/.claude/scripts/worktree-add.sh <branch>` で作る（/tmp に作ると gitleaks の hook が失敗する）
- 実装を始めたあとでテストだけをコミットするときは、テストの部分だけをステージする（`git stash` は他のセッションと
  共有なので使わない）
- 計画と issue だけのブランチは、事実の誤りと矛盾を直したらマージする。言い回しや置き場所の指摘でレビューを
  回し続けない（ADR はマージ後に書き換えられないので、この扱いにしない）

## 方針

- 後方互換を作らない（利用者はまだいない）。形式を変えたら、移行の経路や古い形式の読み込みを残さない。黙った
  フォールバックを足さず、怪しい入力や状態は、推測せずにエラーで止める。既定値は ADR で決めたものだけにする
- NJOY などの参照コードが誤っていそうなときは、誤りを真似ず、正しい方法を実装して `tracing::warn!` を出す。
  理由は ADR に、差の大きさは検証の記録に書く
- 誤差の原因が分かり、おおよそ定量化できたら、記録して詰めるのをやめ、足りない機能を作る
