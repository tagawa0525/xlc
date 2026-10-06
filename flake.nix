{
  description = "xlc: PWR core analysis code (development and verification environment)";

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";

  outputs =
    { nixpkgs, ... }:
    let
      system = "x86_64-linux";
      pkgs = nixpkgs.legacyPackages.${system};
      # 検証用のツール（OpenMC の入力の作成と結果の読み出し、比較のスクリプト）
      python = pkgs.python3.withPackages (
        p: with p; [
          numpy
          h5py
          msgpack
          lxml
          scipy
          pandas
          uncertainties
          setuptools
          setuptools-scm
          pip
        ]
      );
    in
    {
      devShells.${system} = {
        # 開発と検証: Rust のツールと、検証用の Python
        default = pkgs.mkShell {
          packages = with pkgs; [
            cargo
            rustc
            clippy
            rustfmt
            python
            ruff
            poppler-utils
          ];
          # nixpkgs の rustc は標準ライブラリのソースを含まない。rust-analyzer が core と std を
          # 解析できるよう、場所を教える（ないと matches! などが展開できず、誤ったエラーを出す）
          RUST_SRC_PATH = "${pkgs.rustPlatform.rustLibSrc}";
        };
        # 参照用のツール（NJOY, OpenMC）のビルド（docs/verification.md）。
        # OpenMC 0.15.2 が同梱する fmt は GCC 16 でコンパイルできないので GCC 13 を使う
        references = pkgs.mkShell.override { stdenv = pkgs.gcc13Stdenv; } {
          packages = with pkgs; [
            gfortran13
            cmake
            gnumake
            hdf5
            pkg-config
            python
          ];
        };
      };
    };
}
