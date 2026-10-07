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
      rust = pkgs.rust-bin.fromRustupToolchainFile ./rust-toolchain.toml;
      # 検証用の道具（OpenMC の入力の作成と結果の読み出し、比較のスクリプト）
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
        # 開発と検証: Rust の道具と、検証用の Python
        default = pkgs.mkShell {
          packages = with pkgs; [
            rust
            python
            ruff
            poppler-utils
          ];
        };
        # 参照用の道具（NJOY, OpenMC）のビルド（docs/verification.md）。
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
