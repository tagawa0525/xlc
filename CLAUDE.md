# xlc

PWR の解析コード。xs（ENDF から多群定数）→ lattice（MOC の格子計算）→ core（3次元の NEM）の段階からなる。

## 文書

- 決定と理由: docs/adr/（書き方と運用は docs/adr/README.md）
- 段階の目的・範囲・作業の順序・検証の対象: docs/plans/
- 未解決の問題とやること: docs/issues/（書き方と運用は docs/issues/README.md。ブランチ名とコミットの `Refs:`、
  ADR との両方向の対応を含む）
- 検証の記録: examples/*/verification.md。今の結果だけを題目ごとに書き、古い測定や誤りは消して書き直す（経緯は git に
  残る。ADR-0117）。結果ごとに、数値を出した xlc のコミットと日付を書く（「計算」の節。新しい測定には必ず書く。コミットが分からない過去の測定だけ日付）
- コードには How（今の挙動だけ。以前の挙動や変わった経緯は書かず、git、ADR、閉じた issue に任せる）、テストには What、コミットログには Why（決定の理由と採らなかった案は ADR）、コード
  コメントには Why not（その場の、自明に見える別の書き方をしない理由と制約、式の出典と導出の要点。比べた代わりの案は
  ADR か issue に書いて参照する）を書く（ADR-0117）
- 検証の環境（NJOY、OpenMC、核データ）: docs/environment.md。`references/` は git で追跡しない

## 用語と書き方

- 「群」はエネルギー群だけに使う。対称性は「対称性」「対称操作」と書き、コードの識別子でも group を使わない
- 設定（方位角、トラック間隔など）は「今の入力」でなく、それを決めた ADR で示す
- 結論は数値で支えられることだけを書く。確かめていない原因は「候補」「確かめていない」と書く。比べる前に、保存則や
  釣り合い（k と 1/A など）が閉じているか、統計誤差（相関を含む）で区別できる差かを確かめる

## 計算

- 重い計算（OpenMC、NJOY、ライブラリの作り直し、xlc の集合体以上、時間を記録する測定）は、Slurm で 1 本ずつ流す:
  `sbatch --exclusive --wait --export=ALL -J <名前> -o <ログ> --wrap='<コマンド>'`。直接走らせず、`pgrep` で他の計算を見張らない
- OpenMC のスレッド数は docs/environment.md の値を使う
- `examples/<問題>/dvc.yaml` のある例は、`dvc repro examples/<問題>/dvc.yaml`（すべては `dvc repro -P`）で段階を
  走らせ、手で段階を選ばない。重いときは `dvc repro` 全体を上の形で Slurm に流す。走らせたら `dvc.lock` をコミットする（ADR-0123）
- worktree の `references` は main の `references/` への symlink（`.worktreelink`。worktree-add.sh が作る）で、DVC の
  キャッシュと核データを共有する（ADR-0126）
- `examples/vera/vera.msgpack` は git で追跡しない。worktree では main のものをコピーして使い、核種を足して作り直したら
  main の場所にも置く
- 記録の数値には、それを出したコミットと日付を添える（checkout すれば再現できるため）。rebase で hash が変わっても
  測り直さず、道具（remap-commit-refs.sh）で付け替える。文書やテストの変更で hash が最新でなくなっても測り直さない。
  測り直すのは計算（コード、入力、環境）を変えたときだけ（ADR-0122）

## 作業の流れ

- GitHub のリモートはない。PR の代わりに /git-merge でマージする（rebase、チェック、/code-review と指摘への対応、
  `git merge --no-ff`、後片付けを含む）
- チェックは `~/.claude/skills/language-checks/scripts/run-checks.sh` で、終了コードで判定する。`checks | tail && git merge`
  のようにつなぐと、失敗してもマージに進む
- ADR は `design-records check`（`.pre-commit-config.yaml` の hook で、run-checks.sh も `id: design-records` のあるリポジトリで走らせる）でも確かめる
- ADR と issue の番号は `design-records new [--kind issue]` で振る（タグで予約。docs/adr/README.md）。手で最大の番号 + 1 を振らない
- worktree は `~/.claude/scripts/worktree-add.sh <branch>` で作る（/tmp に作ると gitleaks の hook が失敗する）
- 実装を始めたあとでテストだけをコミットするときは、テストの部分だけをステージする（`git stash` は他のセッションと
  共有なので使わない）
- 計画と issue だけのブランチは、事実の誤りと矛盾を直したらマージする。言い回しや置き場所の指摘でレビューを
  回し続けない（ADR はマージ後に書き換えられないので、この扱いにしない）
- 試行は本番のコードに入力ファイルを与えて行う。パッチ、環境変数、入力を組み立てる run shell は作らず、選べない選択肢は
  入力で選べる機能として足す。手元の試行は `trials/`（無視する）に置き、記録に数値を書く試行だけ `git add -f` で
  コミットする。マージの前に残さない試行を消し、残すものは examples/<問題>/ に `git mv` する（ADR-0122）

## 方針

- 後方互換を作らない（利用者はまだいない）。形式を変えたら、移行の経路や古い形式の読み込みを残さない。黙った
  フォールバックを足さず、怪しい入力や状態は、推測せずにエラーで止める。既定値は ADR で決めたものだけにする
- NJOY などの参照コードが誤っていそうなときは、誤りを真似ず、正しい方法を実装して `tracing::warn!` を出す。
  理由は ADR に、差の大きさは検証の記録に書く
- 誤差の原因が分かり、おおよそ定量化できたら、記録して詰めるのをやめ、足りない機能を作る
