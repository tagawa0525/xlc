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
      # ここで `rust` を override しない。flake.nix は DVC のステージの deps でないので、override すると段階をビルドする
      # コンパイラが黙って変わる（#0080）。開発のツールは rust-toolchain.toml にも書かず（Rust の段階の deps）、同じ版の
      # component を横に置く。rust-analyzer は sysroot を `rustc --print sysroot` で見つける（ADR-0135）
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
      # 検証用のツール（OpenMC の入力の作成と結果の読み出し、比較のスクリプト）が使う Python のパッケージ
      pythonPackages =
        p: with p; [
          numpy
          h5py
          msgpack
          lxml
          scipy
          pandas
          uncertainties
        ];
      python = pkgs.python3.withPackages pythonPackages;
      # 計算の段階の駆動（ADR-0123）。pkgs.dvc をそのまま shell に入れると、Python の setup hook が dvc の依存を
      # PYTHONPATH に足し、検証用の Python のパッケージより先に読まれる。dvc の実行ファイルは自分の依存を
      # ラッパーで見つけるので、実行ファイルだけを入れる
      dvc = pkgs.runCommand "dvc-${pkgs.dvc.version}" { } ''
        mkdir -p $out/bin
        ln -s ${pkgs.dvc}/bin/dvc $out/bin/dvc
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
      # zip ごとの ENDF-B-VIII.0_<part>/ を 1 つの木に並べる（references/endf の形。入力の TOML と dvc.yaml が指す）
      endf = pkgs.runCommand "endf-b-viii.0" { nativeBuildInputs = [ pkgs.unzip ]; } ''
        mkdir -p $out
        ${lib.concatMapStrings (zip: "unzip -q ${zip} -d $out\n") (lib.attrValues endfZips)}
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
      referencesPython = pkgs.python3.withPackages (p: pythonPackages p ++ [ openmcPython ]);
    in
    {
      packages.${system} = {
        endf-neutrons = endfZips.neutrons;
        endf-thermal-scatt = endfZips.thermal_scatt;
        endf-decay = endfZips.decay;
        endf-nfy = endfZips.nfy;
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
      };
    };
}
