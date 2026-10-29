# xlc

PWR の解析コード。xs（ENDF から多群定数）→ lattice（MOC の格子計算）→ core（3次元の NEM）のステージからなる。

## 文書

- 決定と理由: docs/adr/（書き方と運用は docs/adr/README.md）
- ステージの目的・範囲・作業の順序・検証の対象: docs/plans/
- 未解決の問題とやること: docs/issues/（書き方と運用は docs/issues/README.md。ブランチ名とコミットの `Refs:`、
  ADR との両方向の対応を含む）
- 検証の記録: examples/*/verification.md（評価済みデータを変えた検証は、混ぜずに examples/<問題>/verification-<略称>.md に分ける。ADR-0175）。今の結果だけを題目ごとに書き、古い測定や誤りは消して書き直す（経緯は git に
  残る。ADR-0117）。結果ごとに、数値を出した xlc のコミットと日付を書く（「計算」の節。新しい測定には必ず書く。コミットが分からない過去の測定だけ日付）
- 検証の記録の表のうち `<!-- table: <名前> -->` 〜 `<!-- /table -->` の間は、比較のステージの JSON から
  `python3 examples/verification_tables.py <記録>` が描く（手で直さない。ADR-0176）。`--check` は hook と `cargo test`
  が実行する。JSON は数十 KB ほどの要約で、DVC の出力（`cache: false`）を git で追跡する（`dvc checkout` なしで走る）
- コードには How（今の挙動だけ。以前の挙動や変わった経緯は書かず、git、ADR、閉じた issue に任せる）、テストには What、コミットログには Why（決定の理由と採らなかった案は ADR）、コード
  コメントには Why not（その場の、自明に見える別の書き方をしない理由と制約、式の出典と導出の要点。比べた代わりの案は
  ADR か issue に書いて参照する）を書く（ADR-0117）
- 検証の環境（NJOY、OpenMC、核データ）: docs/environment.md。`references/` は git で追跡しない

## 用語と書き方

- 「群」はエネルギー群だけに使う。対称性は「対称性」「対称操作」と書き、コードの識別子でも group を使わない
- 設定（方位角、トラック間隔など）は「今の入力」でなく、それを決めた ADR で示す
- 結論は数値で支えられることだけを書く。確かめていない原因は「候補」「確かめていない」と書く。比べる前に、保存則や
  バランス（k と 1/A など）が閉じているか、統計誤差（相関を含む）で区別できる差かを確かめる

## 計算

- 重い計算（OpenMC、NJOY、ライブラリの作り直し、xlc の集合体以上、時間を記録する測定）は、Slurm で 1 本ずつ流す:
  `~/.claude/scripts/slurm-run.sh [--timed] -t <見積もり> <名前> '<コマンド>'`（排他で、終わるまで待つ）。時間を記録する測定は `--timed` を付ける（割り込まれないパーティション timed で、走っている main のジョブを止めて先に走る）。付けないジョブは main で、timed のジョブに止められうるので、その時間は記録に書かない（`TIME:` 行に `preemptible:` の印が付く）。直接の `sbatch` は hook が止める。直接実行せず、`pgrep` で他の計算を監視しない。時間を記録しない長い計算は `-c <CPU 数>`（ハードウェアスレッドで数える）を付けて排他なしに流せるが、その時間は記録に書かない。
- OpenMC のスレッド数は docs/environment.md の値を使う
- lattice のスレッドの数は、`RAYON_NUM_THREADS` が無ければプロセスに許された論理 CPU（SMT を含む）の数で、r995 の排他のジョブでは
  32（1 コア 2 スレッドの 16 コア。`lscpu`）。始めにログに出す（`thread pool threads=`）。xs の `--threads` の既定は
  機械の物理コアの数をプロセスに許された CPU の数で抑えたもの、runner の予算は許された CPU の中の物理コアの数で、r995 の排他のジョブでは
  どちらも 16、`-c 16` では xs 16、runner 8（ADR-0214、ADR-0193）。xs の数はログに出ないので、記録と
  依頼に書くときは、この規則と、ジョブで付けた `--threads`、`-c` の数から書く（#0186）
- 数時間の OpenMC の参照解は AWS のスポット（`aws_run.sh --type auto`。ADR-0171）、xlc の計算は条件が変わりやすいので r995 の
  Slurm、短い OpenMC はローカル。機種の選び方と、オンデマンドと m7a 系を使わないことは docs/environment.md の「外部の計算機（AWS）」（ユーザー 2026-10-13、2026-10-14）
- 重い計算を足す前に、既存の結果で足りるか、結論が変わりうるかを見る。形式（履歴、テーブル、ライブラリ）を変える変更の前に、
  読めなくなる出力と作り直す時間を数えて依頼と報告に書く。読めなくなった出力の再実行は、それを読む作業が次に来るまで後回し（ユーザー 2026-10-11、2026-10-24）
- 高速化の途中は、高速化の測定以外のジョブを流さない。遅い bin で流した計算は流し直しになる。後回しは 0 コスト、見込みの外れた
  ジョブは理由を書いて止める（ユーザー 2026-10-24）
- 仕様書に KENO-VI の参照解がある問題は ENDF/B-VII.0 の検証用ライブラリで KENO と直接比べ、OpenMC の参照解は仕様書に参照のない
  問題（燃焼、分岐、HFP）と診断だけに作る（ADR-0175）。3D の全炉心と Problem 4 の OpenMC は作らない（ADR-0124、docs/plans/m3-core-hzp.md。ユーザー 2026-10-15）
- 結果を変えない deps の変更で、重いステージを再実行せず `dvc commit` で付け直してよいのは、ADR-0212 の (b) の条件（値の同じ TOML だけ）を
  満たすときと、ADR が決めた移行（ADR-0159、ADR-0210）だけ。参照のステージ（NJOY、OpenMC）が使わない Python のパッケージは、flake.nix の `pythonPackages` にだけ足す（`references` の
  devShell の store のパスが参照解の deps。食い違いは pre-commit の hook `references-shell` が止める。ADR-0212）
- `examples/<問題>/dvc.yaml` のある例は、`dvc repro examples/<問題>/dvc.yaml`（すべては `dvc repro -P`）でステージを
  実行する。重いときは全体を上の形で Slurm に流す。重いステージを 1 本ずつ流す、または参照解の待ちで全体が失敗するときは、
  目標を指定する `dvc repro <yaml>:<ステージ>`（上流の changed は実行される。ADR-0184）。`-s` は上流を飛ばすので、rebase のあとに
  数秒のステージを再実行する /dvc-relock だけに使う。実行したら `dvc.lock` をコミットする（上流の `build` と `python-ast` が書き換えたリポジトリの直下の `dvc.lock` と `python-ast/` も。ADR-0159、ADR-0210）
- 核種ごとの xs のステージ（`xs-nuclide`）は、`dvc repro` でなく `python3 examples/stage_runner.py <yaml> xs-nuclide` で流す（changed の
  ものを CPU 時間の長い順に同時に実行して `dvc commit` する。上流の changed も実行する。下流の `xs` は続けて `dvc repro <yaml>:xs`。ADR-0193）。
  runner が書き直す CPU 時間の記録 `<yaml のディレクトリ>/stage-cpu-seconds.json` も `dvc.lock` と一緒にコミットする
- `dvc repro` を投げる前に、走りうるステージを数えて見積もる（出力を戻した、ジョブを実行する worktree で）。
  changed の数は `dvc status` の文字列を grep で数えず（タブで始まる `changed deps:` の行を拾って誤り、2026-10-28 に 468 を 942 と数えた）、
  `scripts/count-changed-stages.py [<yaml>...]` で数える。
  `dvc status --with-deps <yaml>:<ステージ>` の changed と、目標の上流（`dvc dag --dot <yaml>:<ステージ>`）のうちその下流に
  あるもの（`dvc status` と `dvc repro --dry` は、上流が出力を書き換えた後に走る下流を出さない。上流が同じ出力を書けば下流は
  飛ぶので上限）。走るステージに参照解（数時間の OpenMC）が入っていれば、投げずに原因（deps の文書だけの変更など）を調べる。
  時間は、過去のログで実行したステージの数の一覧（`grep -o "Running stage '[^'@]*" <ログ> | sort | uniq -c`）の全体が同じジョブから取る。
  2026-10-22 に、65 核種の作り直しを多群の作成のないジョブから 45 分と見積もって 2 時間 5 分かかり（Slurm 452。ラダーの和をベクトル化した ADR-0188 の後は 38 分、Slurm 473）、文書だけを
  直した deps で changed になった 5 時間の OpenMC を目標の上流に入れて投げ、12 分で止めた（Slurm 454）
- worktree の `references` は main の `references/` への symlink（`.worktreelink`。worktree-add.sh が作る）で、DVC の
  キャッシュと核データを共有する（ADR-0126）
- `examples/vera/vera.msgpack` と VERA の参照解（OpenMC のライブラリ、statepoint、k）は examples/vera/dvc.yaml のステージの
  出力で、git で追跡しない。worktree では `scripts/worktree-restore-outputs.sh`（中身は `dvc checkout`）でキャッシュから戻す（ADR-0130）。ステージが呼ぶ
  `bin/` の実行ファイルは、リポジトリの直下の `dvc.yaml` の `build` の出力で、同じスクリプトで戻る（ADR-0159）。
  ステージの deps の .py は、元のファイルでなく、直下の `dvc.yaml` の `python-ast` が書く正規形（docstring とコメントを除いた形。
  `python-ast/<パス>.ast`、git で追跡）にする（cmd は元の .py を呼ぶ。ADR-0210）。.py を deps に足すときは、`python-ast` の
  foreach にも足す（`scripts/test_normalize_python.py` が食い違いを止める）。
  NJOY と OpenMC を使うステージは cmd が `nix develop ../..#references -c` で包んであるので、`dvc repro` は default の
  devShell で実行する
- 記録の数値には、それを出したコミットと日付を添える（checkout すれば再現できるため）。ブランチの上で実行した計算は、実行した
  コミットと分岐点（main のコミット）を書く。rebase で deps が変わっても再実行せず、分岐点から復元する。数秒で終わる
  ステージだけは再実行してよい（ADR-0130）。rebase で hash が変わっても
  測り直さず、ツール（`commit-refs remap`。rebase と amend では、`remap.commitRefs` が true なら post-rewrite hook が自動で呼ぶ（失敗したときだけ手で呼ぶ）。付け替えはコミットされないので、自分でコミットする）で付け替える。文書やテストの変更で hash が最新でなくなっても測り直さない。
  再測定するのは計算（コード、入力、環境）を変えたときだけ（ADR-0122）

## 作業の流れ

- GitHub のリモートはない。PR の代わりに /git-merge でマージする（rebase、チェック、/code-review と指摘への対応、
  `git merge --no-ff`、後片付けを含む）
- `dvc.yaml` のある例を含むブランチは、/git-merge の rebase のあとに /dvc-relock、マージのあとに同スキルの後半
- 1 つのトピックの始めからマージまでは /topic-start に従う。委任を受けた担当も作業の前に読む（記録用の長い計算の前に
  rebase とレビューを済ませるなど、/delegate の手順の順と違うところがある）
- チェックは `~/.claude/skills/language-checks/scripts/run-checks.sh` で、終了コードで判定する。`checks | tail && git merge`
  のようにつなぐと、失敗してもマージに進む
- 差分が文書（`.md`）だけのとき、run-checks.sh は fmt、clippy、全体の `cargo test` の代わりに docs クレートのテスト
  （`cargo test -p xlc-docs`）だけを実行する。宣言は `scripts/run-rust-docs-tests.sh`（nixfiles の ADR-0057）。ほかの差分では
  全体の後にも宣言を実行し、宣言に実行権限が無ければ差分によらず失敗する
- Python のテスト（追跡している `test_*.py` と `*_test.py`）は `scripts/run-python-tests.sh` が python-tests の devShell（references の Python に検証用のツールのパッケージを足したもの。ADR-0212）で実行し、run-checks.sh が
  それを呼ぶ（nixfiles の ADR-0036）。run-checks.sh は Python のステージを .py の変更があるときだけ実行するので、
  手元では `scripts/run-python-tests.sh` で全部を実行できる
- ADR は `design-records check`（`.pre-commit-config.yaml` の hook で、run-checks.sh も `id: design-records` のあるリポジトリで実行する）でも確かめる
- ADR と issue の番号は `design-records new [--kind issue]` で振る（タグで予約。docs/adr/README.md）。手で最大の番号 + 1 を振らない
- worktree は `~/.claude/scripts/worktree-add.sh <branch>` で作る（/tmp に作ると gitleaks の hook が失敗する）
- 実装を始めたあとでテストだけをコミットするときは、テストの部分だけをステージする（`git stash` は他のセッションと
  共有なので使わない）
- マージしないで残すブランチ（ユーザーの判断で消さない）は docs/branches.md
- 計画と issue だけのブランチは、事実の誤りと矛盾を直したらマージする。言い回しや置き場所の指摘でレビューを
  回し続けない（ADR の凍結する部分は、マージ後に事実の誤りを直せないので、この扱いにしない。docs/adr/README.md の「運用」）
- 試行は本番のコードに入力ファイルを与えて行う。パッチ、環境変数、入力を組み立てる run shell は作らず、選べない選択肢は
  入力で選べる機能として足す。手元の試行は `trials/`（無視する）に置き、記録に数値を書く試行だけ `git add -f` で
  コミットする。マージの前に残さない試行を消し、残すものは examples/<問題>/ に `git mv` する（ADR-0122）

## 方針

- 後方互換を作らない（利用者はまだいない）。形式を変えたら、移行の経路や古い形式の読み込みを残さない。黙った
  フォールバックを足さず、怪しい入力や状態は、推測せずにエラーで止める。既定値は ADR で決めたものだけにする
- NJOY などの参照コードが誤っていそうなときは、誤りを真似ず、正しい方法を実装して `tracing::warn!` を出す。
  理由は ADR に、差の大きさは検証の記録に書く
- 誤差の原因が分かり、おおよそ定量化できたら、記録して詰めるのをやめ、足りない機能を作る
- 誤差の原因の候補や次に詰める差を並べる前に（宛先によらない。依頼、ADR の検討した案、ユーザーへの提案も同じ）、その量の語
  （移行、散乱行列、群の重みなど）で open の issue と ADR の「候補」「確かめていない」「後回し」の行を検索し、既に挙がっている
  候補を一覧に入れる。差を挙げるときは、検証の記録の値と ± で統計誤差と区別できるかを先に確かめ、issue にある見方を添える。
  2026-10-27 に、#0005、ADR-0076、ADR-0133 に残っていた散乱行列の群の中の重み（1/E）を ADR-0202 の候補から落とし、
  散乱核を ADR-0203〜0205 で先に詰めた（原因は ADR-0207 で分かった）。2026-10-28 に、1C の 0 の点の差の 1 リングと 4 リングの
  違い（+53 ± 22 と +102 ± 20 pcm）をユーザーへの次の手の 1 番目に挙げ、比較の JSON で xlc どうしの差が +17 pcm、残りは
  OpenMC の統計の範囲（#0133 の経緯も揺らぎの大きさの中と書いていた）と分かって外した
- 判断の軸は正確さ > 実行速度 > 保守の容易さ。自己遮蔽は UFG のまま、物理的に正しい改善で速くし、衝突確率の精度（1e-9）は保つ（ユーザー 2026-10-12）。ここでの「速度」は実行速度で、ライブラリの作り直しや再実行を避ける意味ではない（無駄に CPU を使わないことは、重い計算の見積もりと後回しの規則が別に扱う）。保守性を捨てているのではなく、3 つが衝突したときの順で、保守性が上がる変更（複製を共通のクレートにまとめるなど）は、実行ファイルがどのみち変わるときなどに費用を数えて行う（ユーザー 2026-10-28）
- 速度のための近似は入力の項目として足し、既定は物理的に正しい方法にする（参照が数時間で既定にできなければ、参照からの差を測った
  速い設定。ADR-0209）。厳密な経路は入力で選べる形で必ず残す（ほかの手法の参照解にする。ユーザー 2026-10-27）。本番（記録に載せる燃焼、
  p2 のステージ）は既定、検討（結論が ±20 pcm で変わらないもの）は `shielding = "step-start"` などの簡易でよい（ADR-0197、ユーザー 2026-10-24）
- 問題の大きさに釣り合わない時間は性能の欠陥として issue にし、同じ計算を繰り返す前に最小の入力で原因を測って直す。計画の見積もりにも
  当てはめる。商用コードの桁から数倍までは手法の代価として許容、2 桁は設計の誤り（ユーザー 2026-10-13、2026-10-24。例: 1 次元の平板の 55 分 #0069、
  1 集合体のテーブル 110 時間、ADR-0209 の後の自己遮蔽 #0180）
- 燃焼 → 炉心の燃焼 → 動特性を最小の入力（ピンセル、集合体、小さい炉心）で縦に通して時間の内訳を測り、ネックを掴んでから詰める。
  大きい計算（9 種のテーブル、全炉心）は精度の原因を詰めて直した後。簡単に直せる非効率はすぐ直す（ユーザー 2026-10-14、2026-10-24）
- 環境整備（DVC、チェック、スキル、references の片付け、リファクタリング）も issue に書いて計画的に進める（ユーザー 2026-10-12）
