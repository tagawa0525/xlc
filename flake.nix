{
  description = "xlc: PWR core analysis code (development and verification environment)";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    # nixpkgs の Rust は stable の公開から遅れるので、rust-toolchain.toml の版を rust-overlay から取る
    rust-overlay = {
      url = "github:oxalica/rust-overlay";
      inputs.nixpkgs.follows = "nixpkgs";
    };
  };

  outputs =
    { nixpkgs, rust-overlay, ... }:
    let
      system = "x86_64-linux";
      pkgs = import nixpkgs {
        inherit system;
        overlays = [ rust-overlay.overlays.default ];
      };
      inherit (pkgs) lib;
      rust = pkgs.rust-bin.fromRustupToolchainFile ./rust-toolchain.toml;
      # ここで `rust` を override しない。override するとステージの実行ファイルをビルドするコンパイラが変わる（#0080）。
      # 開発のツールは rust-toolchain.toml にも書かず（DVC の build のステージの deps）、同じ版の component を横に置く。
      # rust-analyzer は sysroot を `rustc --print sysroot` で見つける（ADR-0135）
      rustAnalyzerBin = rust.availableComponents.rust-analyzer;
      # rust-analyzer の flycheck（cargo check）と build script の成果物を target/rust-analyzer に分ける（ADR-0147）。
      # `cargo.targetDir` は rust-analyzer.toml では効かず（#0079）、クライアントの設定は git で追跡できないので、
      # 起動する rust-analyzer の環境で CARGO_TARGET_DIR を決める。run-checks.sh や cargo を直接実行する shell には
      # 設定しない。git の外で起動したら git が失敗して止まる
      rustAnalyzer = pkgs.writeShellScriptBin "rust-analyzer" ''
        set -euo pipefail
        CARGO_TARGET_DIR="$(${pkgs.git}/bin/git rev-parse --show-toplevel)/target/rust-analyzer"
        export CARGO_TARGET_DIR
        exec ${rustAnalyzerBin}/bin/rust-analyzer "$@"
      '';
      # 参照のステージ（references の devShell で走る OpenMC と NJOY のスクリプト）が使う Python のパッケージ。
      # references の devShell の store のパスはそのステージの deps なので、ここを変えると参照解が changed になる。
      # 参照のステージが使わないパッケージは pythonPackages にだけ足す（#0183）
      referencesPythonPackages =
        p: with p; [
          numpy
          h5py
          msgpack
          lxml
          scipy
          pandas
          uncertainties
        ];
      # 検証用のツール（比較のスクリプト、dvc.lock を読むスクリプト）とテストが使う Python のパッケージ
      pythonPackages = p: referencesPythonPackages p ++ [ p.pyyaml ];
      python = pkgs.python3.withPackages pythonPackages;
      # 計算のステージの駆動（ADR-0159）。pkgs.dvc をそのまま shell に入れると、Python の setup hook が dvc の依存を
      # PYTHONPATH に足し、検証用の Python のパッケージより先に読まれる。pkgs.dvc の bin/dvc（makeWrapper）も
      # 使わない: PATH の先頭に素の python3 を足し、ステージの cmd の python3 が検証用の Python でなくなる（#0101）。
      # 中の .dvc-wrapped は自分の依存の site のディレクトリを自分で足すので、PATH を変えずに直接呼ぶ。
      # .dvc-wrapped は nixpkgs の Python のラッパーの中の名前なので、無ければビルドで止める
      #
      # DVC は vars と params の TOML を tomlkit で読む。nixpkgs の tomlkit 0.15.0 は、共通の接頭辞を持つ表の見出しで読む時間が
      # 見出しの数の 2 乗で延びる（tomlkit の issue 479。0.15.1 で線形になった）ので、DVC の Python だけ直った版にする（#0268）。
      # nixpkgs の tomlkit が 0.15.1 以上になったら、この override は評価で止まるので消す
      dvcTomlkit =
        assert lib.assertMsg (lib.versionOlder pkgs.python3Packages.tomlkit.version "0.15.1")
          "nixpkgs の tomlkit が 0.15.1 以上になった。flake.nix の dvcTomlkit の override を消す（#0268）";
        pkgs.python3Packages.tomlkit.overridePythonAttrs rec {
          version = "0.15.2";
          src = pkgs.fetchPypi {
            pname = "tomlkit";
            inherit version;
            hash = "sha256-g+MnVbodINk0C9Bp7pfp9iMrojRrY//LpCXeAffjJ6Y=";
          };
        };
      dvcPackage = pkgs.python3Packages.toPythonApplication (
        pkgs.python3Packages.dvc.override { tomlkit = dvcTomlkit; }
      );
      dvc = pkgs.runCommand "dvc-${dvcPackage.version}" { } ''
        test -x ${dvcPackage}/bin/.dvc-wrapped || {
          echo "dvc のパッケージ（dvcPackage）に bin/.dvc-wrapped が無い（#0101）" >&2
          exit 1
        }
        mkdir -p $out/bin
        cat > $out/bin/dvc <<EOF
        #!${pkgs.runtimeShell}
        export PYTHONNOUSERSITE=true
        exec ${dvcPackage}/bin/.dvc-wrapped "\$@"
        EOF
        chmod +x $out/bin/dvc
      '';

      # 検証用の資産（ADR-0124）。外から取るものは URL と sha256 で固定する

      # 評価済み核データ ENDF/B-VIII.0。sha256 は NNDC の zip そのもののハッシュで、`sha256sum` で照らせるように
      # fetchzip（展開した木の NAR ハッシュ）でなく fetchurl で取る
      endfZip =
        part: sha256:
        pkgs.fetchurl {
          url = "https://www.nndc.bnl.gov/endf-b8.0/zips/ENDF-B-VIII.0_${part}.zip";
          inherit sha256;
        };
      endfZips = {
        neutrons = endfZip "neutrons" "a1cced11710eb36690c5ad7b64d414465c5849f5428904049c864e9beca28894";
        thermal_scatt = endfZip "thermal_scatt" "042552d380e7fc3a1e7db948523b317ebaff260a08f5c8a786e94918e6e838f8";
        decay = endfZip "decay" "9087395f6d74fe35b47bb4c8b37ef53b44e67569196c933f835a53cdf8d03d5e";
        nfy = endfZip "nfy" "92c5371fdb21eecf4989f48828671b904186abc6386b3d7510c8fcee2ee5ffcf";
      };
      # 検証用の ENDF/B-VII.0（仕様書の KENO-VI と核データの版を揃えて比べる。ADR-0175）と、版ごとの燃焼チェーンの decay と
      # nfy（#0303）。NNDC の zip はディレクトリを持たない平らな形なので、VIII.0 と同じ ENDF-B-VII.0_<part>/ に展開する
      endfViiZip =
        part: sha256:
        pkgs.fetchurl {
          url = "https://www.nndc.bnl.gov/endf-b7.0/zips/${part}-ENDF-VII0.endf.zip";
          inherit sha256;
        };
      endfViiZips = {
        neutrons = endfViiZip "n" "935a1c8156958072df7363239f8c7a6a514e1d20a4053986549b9b4e093f23b9";
        thermal_scatt = endfViiZip "tsl" "754d4e41280a7339b9c352378c7b16dd659d9dad244dd4c68dafa444e86f160e";
        decay = endfViiZip "dec" "0778766165d723ddfa41c44a3f3bdb46d15c631c57364196d7cf012a65b84046";
        nfy = endfViiZip "nfy" "13f5589462965a572c582af5129ff9c903e12e4f43bdc6b156a958296cdbaf96";
      };
      # 版ごとの燃焼チェーンとライブラリ（#0303）の ENDF/B-VII.1。zip は <part>/ のディレクトリを持つ
      endfVii1Zip =
        part: sha256:
        pkgs.fetchurl {
          url = "https://www.nndc.bnl.gov/endf-b7.1/zips/ENDF-B-VII.1-${part}.zip";
          inherit sha256;
        };
      endfVii1Zips = {
        neutrons = endfVii1Zip "neutrons" "ed106bfddb537cf46ac944536baec3b03b08d92a8efd97ff4e782531a82f99f3";
        thermal_scatt = endfVii1Zip "thermal_scatt" "796629f66c4a8f674fc1b1040c2336ce42152ed324d93d6be60d81e4bca1aff2";
        decay = endfVii1Zip "decay" "25edd580cd36533194d43ae71503b7e2eebe29141f973c1621e2b026bb671df8";
        nfy = endfVii1Zip "nfy" "6338f5466ede3fee50008adc3d24af0923df6ac051800be133c34a210c0f54c3";
      };
      # 版ごとの燃焼チェーンとライブラリ（#0303）の ENDF/B-VIII.1。tar.gz は <part>-version.VIII.1/ のディレクトリを持つ
      endfViii1Tar =
        part: sha256:
        pkgs.fetchurl {
          url = "https://www.nndc.bnl.gov/endf-releases/releases/B-VIII.1/${part}/${part}-version.VIII.1.tar.gz";
          inherit sha256;
        };
      endfViii1Tars = {
        neutrons = endfViii1Tar "neutrons" "decff90016bfb5c25c8d4c7fb8d81f94ff2f104853165f0fe6a21bf61c4164e4";
        thermal_scatt = endfViii1Tar "thermal_scatt" "bfb824a9719253e4be6e9be5b1faf76b1b78452b941005babddadb3294c9a184";
        decay = endfViii1Tar "decay" "8639581a4882b35855f7e98e14e537470c6b4b2a6d0673b8ac982baaa50a9d79";
        nfy = endfViii1Tar "nfy" "2cf9b251c66fb0d8da96d03ad017fc3893cd3387a1d5f8cec37223c76319716a";
      };
      # 版と部分ごとの ENDF-B-<版>_<part>/ を 1 つの木に並べる（references/endf の形。入力の TOML と dvc.yaml が指す）。
      # 中のディレクトリの名前が版で違うので、VII.1 と VIII.1 は一時の場所に展開してから名前を揃える
      endf = pkgs.runCommand "endf" { nativeBuildInputs = [ pkgs.unzip ]; } ''
        mkdir -p $out
        ${lib.concatMapStrings (zip: "unzip -q ${zip} -d $out\n") (lib.attrValues endfZips)}
        ${lib.concatStrings (
          lib.mapAttrsToList (part: zip: "unzip -q ${zip} -d $out/ENDF-B-VII.0_${part}\n") endfViiZips
        )}
        ${lib.concatStrings (
          lib.mapAttrsToList (
            part: zip: "unzip -q ${zip} -d unpack-vii1 && mv unpack-vii1/${part} $out/ENDF-B-VII.1_${part}\n"
          ) endfVii1Zips
        )}
        ${lib.concatStrings (
          lib.mapAttrsToList (part: tar: ''
            mkdir -p unpack-viii1 && tar -xzf ${tar} -C unpack-viii1
            mv unpack-viii1/${part}-version.VIII.1 $out/ENDF-B-VIII.1_${part}
          '') endfViii1Tars
        )}
      '';

      # NJOY2016 2016.79 に、ACER の LAW=0 の書き出しを直すパッチ（docs/environment.md）を当てる
      njoy2016 = pkgs.gcc13Stdenv.mkDerivation {
        pname = "njoy2016";
        version = "2016.79";
        src = pkgs.fetchFromGitHub {
          owner = "njoy";
          repo = "NJOY2016";
          # タグ 2016.79
          rev = "ac5adf5f33d893e42f2eed7fb286b0d51c7580da";
          hash = "sha256-2eoikGtQZU85n6mue0cl0YkCaLDuGDrpmat8H0d7JJo=";
        };
        patches = [ ./docs/patches/njoy2016-acer-law0.patch ];
        nativeBuildInputs = with pkgs; [
          cmake
          gfortran13
          python3
        ];
        # CMakeLists.txt に install の規則がない
        installPhase = ''
          runHook preInstall
          install -Dm755 njoy $out/bin/njoy
          runHook postInstall
        '';
      };

      openmcVersion = "0.15.2";
      openmcRev = "e23760b0264c66fb7bb373aa0596801e5209d920"; # タグ v0.15.2
      # 同梱の fmt、pugixml、xtensor をサブモジュールから使う（手でビルドしたときと同じ）
      openmcSrc = pkgs.fetchFromGitHub {
        owner = "openmc-dev";
        repo = "openmc";
        rev = openmcRev;
        fetchSubmodules = true;
        hash = "sha256-BsgrN8pjNwIvaRVHjAMYNns1IfJ6He876JrEKnwd+MU=";
      };
      # OpenMC の CMake は版を git describe か .git_archival.txt から読む。取ったソースには .git がなく、
      # .git_archival.txt は git archive が埋めるテンプレートのままなので、埋めた形に書き換える
      openmcArchival = ''
        cat > .git_archival.txt <<ARCHIVAL
        commit: ${openmcRev}
        describe-name: v${openmcVersion}
        ARCHIVAL
      '';
      # 同梱の fmt が GCC 16 でコンパイルできないので GCC 13 を使う
      openmc = pkgs.gcc13Stdenv.mkDerivation {
        pname = "openmc";
        version = openmcVersion;
        src = openmcSrc;
        postPatch = openmcArchival;
        nativeBuildInputs = with pkgs; [
          cmake
          pkg-config
        ];
        buildInputs = [ pkgs.hdf5 ];
        cmakeFlags = [
          "-DOPENMC_USE_OPENMP=on"
          # テストの実行ファイルは使わない
          "-DOPENMC_BUILD_TESTS=off"
        ];
      };

      # OpenMC の Python API（openmc.data で NJOY から HDF5 を作る、モデルを書く、statepoint を読む）と、
      # それが使う endf（PyPI の sdist から作る）
      endfPython = pkgs.python3.pkgs.buildPythonPackage rec {
        pname = "endf";
        version = "0.1.12";
        pyproject = true;
        src = pkgs.fetchPypi {
          inherit pname version;
          hash = "sha256-l8obeiSPzCbFip9qgl5ItmQ8im4RFzNECtoK5CXUWmg=";
        };
        build-system = with pkgs.python3.pkgs; [
          setuptools
          setuptools-scm
          pybind11
        ];
        dependencies = with pkgs.python3.pkgs; [
          numpy
          uncertainties
        ];
        pythonImportsCheck = [ "endf" ];
      };
      openmcPython = pkgs.python3.pkgs.buildPythonPackage {
        pname = "openmc";
        version = openmcVersion;
        pyproject = true;
        src = openmcSrc;
        # openmc.lib は、パッケージの中の libopenmc.so を読み込む。symlink にしても wheel を作るときに
        # コピーされるので、コピーして置く
        postPatch = ''
          cp ${openmc}/lib/libopenmc.so openmc/lib/
        '';
        # setuptools-scm は .git がないと版を決められない
        env.SETUPTOOLS_SCM_PRETEND_VERSION = openmcVersion;
        build-system = with pkgs.python3.pkgs; [
          setuptools
          setuptools-scm
        ];
        dependencies =
          with pkgs.python3.pkgs;
          [
            numpy
            h5py
            scipy
            pandas
            lxml
            uncertainties
            setuptools
          ]
          ++ [ endfPython ];
        # 対話と作図にだけ使う依存で、検証のスクリプトは使わない
        pythonRemoveDeps = [
          "ipython"
          "matplotlib"
        ];
        pythonImportsCheck = [ "openmc" ];
      };
      referencesPython = pkgs.python3.withPackages (p: referencesPythonPackages p ++ [ openmcPython ]);
      # Python のテスト（scripts/run-python-tests.sh）: 検証用のツールのパッケージと OpenMC の Python API の両方を使う
      testsPython = pkgs.python3.withPackages (p: pythonPackages p ++ [ openmcPython ]);
    in
    {
      checks.${system} = {
        # ステージの cmd の python3 が、dvc の PATH でなく検証用の Python になること（#0101）
        dvc-stage-python = pkgs.runCommand "dvc-stage-python" { } ''
          export PATH=${dvc}/bin:${python}/bin:${pkgs.coreutils}/bin HOME=$TMPDIR DVC_SITE_CACHE_DIR=$TMPDIR/dvc-site
          export DVC_NO_ANALYTICS=1
          cd $TMPDIR
          dvc init --no-scm
          dvc stage add -n t -o out.txt "python3 -c 'import msgpack, numpy' && echo ok > out.txt"
          dvc repro
          cp out.txt $out
        '';
        # 自己完結な実行の VM の user-data（examples/aws/aws_user_data.sh。ADR-0155）
        aws-user-data = pkgs.runCommand "aws-user-data" { } ''
          export PATH=${pkgs.bash}/bin:${pkgs.coreutils}/bin:${pkgs.gnugrep}/bin:${pkgs.gnused}/bin
          bash ${./examples/aws/test_aws_user_data.sh} ${./examples/aws/aws_user_data.sh} > $out
        '';
        # スポットの vCPU の残りの読み取りと、1 本の上限（偽の aws。ADR-0171）
        aws-vcpu = pkgs.runCommand "aws-vcpu" { } ''
          export PATH=${pkgs.bash}/bin:${pkgs.coreutils}/bin:${pkgs.gnugrep}/bin:${pkgs.gnused}/bin:${pkgs.gawk}/bin
          bash ${./examples/aws/test_aws_vcpu.sh} ${./examples/aws/aws_common.sh} > $out
        '';
        # 監視 aws_watch.sh の二重起動の防止（flock）と、欠けたジョブをまとめた起動し直し（偽の aws_run.sh。ADR-0171）
        aws-watch = pkgs.runCommand "aws-watch" { } ''
          export PATH=${pkgs.bash}/bin:${pkgs.coreutils}/bin:${pkgs.gnugrep}/bin:${pkgs.gnused}/bin:${pkgs.util-linux}/bin
          bash ${./examples/aws/test_aws_watch.sh} ${./examples/aws/aws_watch.sh} > $out
        '';
        # 実行の定義（examples/<問題>/aws.toml。ADR-0174）を読む部品と、リポジトリの定義
        aws-definition = pkgs.runCommand "aws-definition" { } ''
          export PATH=${pkgs.bash}/bin:${pkgs.coreutils}/bin:${pkgs.gnugrep}/bin:${pkgs.gnused}/bin
          bash ${./examples/aws/test_aws_definition.sh} ${./examples/aws/aws_common.sh} \
            ${./examples/vera/aws.toml} ${./examples/depletion/aws.toml} ${./examples/depletion/p2/aws.toml} ${./examples/aws/aws-check/aws.toml} > $out
        '';
      };

      packages.${system} = {
        endf-neutrons = endfZips.neutrons;
        endf-thermal-scatt = endfZips.thermal_scatt;
        endf-decay = endfZips.decay;
        endf-nfy = endfZips.nfy;
        endfb70-neutrons = endfViiZips.neutrons;
        endfb70-thermal-scatt = endfViiZips.thermal_scatt;
        endfb70-decay = endfViiZips.decay;
        endfb70-nfy = endfViiZips.nfy;
        endfb71-neutrons = endfVii1Zips.neutrons;
        endfb71-thermal-scatt = endfVii1Zips.thermal_scatt;
        endfb71-decay = endfVii1Zips.decay;
        endfb71-nfy = endfVii1Zips.nfy;
        endfb81-neutrons = endfViii1Tars.neutrons;
        endfb81-thermal-scatt = endfViii1Tars.thermal_scatt;
        endfb81-decay = endfViii1Tars.decay;
        endfb81-nfy = endfViii1Tars.nfy;
        inherit endf njoy2016 openmc;
        openmc-python = referencesPython;
      };

      devShells.${system} = {
        # 開発と検証: Rust のツールと、検証用の Python
        default = pkgs.mkShell {
          packages = [
            rust
            rustAnalyzer
            python
            dvc
            pkgs.ruff
            pkgs.poppler-utils
          ];
          # stdenv が RUNPATH に足す、devShell に入ったディレクトリの outputs/out/lib を止める。足すと同じソースでも
          # worktree ごとに実行ファイルのバイトが変わり、DVC のステージが再実行する（ADR-0159）。そこに置くライブラリはない
          NIX_NO_SELF_RPATH = "1";
          # 外の shell の PYTHONPATH を引き継ぐと、検証用の Python と DVC のステージの python3 が、そのパッケージを site-packages
          # より先に読み、検証用でない版で黙って走る（#0103）。この devShell の Python は withPackages で、PYTHONPATH を要らない
          shellHook = ''
            unset PYTHONPATH
          '';
        };
        # 参照の計算: NJOY、OpenMC と、OpenMC の Python API の入った Python（docs/environment.md）
        references = pkgs.mkShell {
          packages = [
            njoy2016
            openmc
            referencesPython
          ];
          # 外の shell（direnv で入った default の devShell など）の PYTHONPATH を引き継ぐと、その shell の Python の
          # パッケージが referencesPython のものより先に読まれる。referencesPython は PYTHONPATH なしで自分の
          # パッケージを見つける
          shellHook = ''
            unset PYTHONPATH
          '';
        };
        # Python のテスト: references の Python に検証用のツールのパッケージを足したもの。参照のステージの deps に
        # 入らないので、検証用のパッケージを足しても参照解は changed にならない
        python-tests = pkgs.mkShell {
          packages = [ testsPython ];
          shellHook = ''
            unset PYTHONPATH
          '';
        };
      };
    };
}
