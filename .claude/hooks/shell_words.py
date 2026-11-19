"""Claude Code の PreToolUse の hook が共有する、bash のコマンドの文字列の読み方（ADR-0223、ADR-0257）。

コマンドの文字列を実行される単純なコマンドの語の並びに分け（`simple_commands`）、先頭の代入と包むコマンドを外し
（`unwrap`）、dvc のサブコマンドを直接呼ぶかを見る（`runs_dvc`）。読み分ける文法と読まないものは ADR-0223 の決定と帰結。

hook はシェルから devShell の外で実行されるので、標準ライブラリだけを使う。
"""

import os
import re

# 単純なコマンドを区切る記号。改行（ヒアドキュメントの本文を飛ばす）とリダイレクト（< と >）は別に扱う
SEPARATORS = set(";&|")
ASSIGNMENT = re.compile(r"[A-Za-z_][A-Za-z0-9_]*=")
SHELLS = {"bash", "sh", "zsh", "dash"}
# 次の語をコマンドとして実行する予約語（`for`、`case` などの後の語はコマンドでない）
RESERVED = {"if", "then", "elif", "else", "while", "until", "do", "{", "!"}
PYTHONS = re.compile(r"python(3(\.\d+)?)?")


class _Parser:
    """bash のコマンドの文字列を、実行される単純なコマンドの語の並びに分ける。

    引用（`'…'`、`"…"`、`$'…'`、`\\`）、コメント、ヒアドキュメントの本文（コマンドとして実行されない。区切りの語を
    引用しない本文の `$(…)` とバッククォートは bash が実行するが、読まない）を読み分け、`$(…)`、`<(…)`、
    `>(…)` と `"…"` の中の `$(…)` は中のコマンドも並びに入れる。バッククォート、算術式、`case` の `)` は読み分けない（区切りとして扱う）。
    """

    def __init__(self, text):
        self.s = text
        self.i = 0
        self.commands = []
        self.heredocs = []  # (区切りの語, <<- か)

    def parse(self, close_paren=False):
        """close_paren なら対応する `)` で戻る。それ以外は文字列の終わりまで読む。"""
        words = []
        word = None  # None は語が始まっていない
        redirect_target = False
        depth = 0

        def end_word():
            nonlocal word, redirect_target
            if word is not None:
                if redirect_target:
                    redirect_target = False
                else:
                    words.append(word)
            word = None

        def end_command():
            # コマンドの区切りと置換の始まりで、リダイレクトの先の語を待つのをやめる（`<(…)`、`> $(…) && …` の後の
            # 最初の語はコマンド名）
            nonlocal words, redirect_target
            end_word()
            redirect_target = False
            if words:
                self.commands.append(words)
            words = []

        s = self.s
        while self.i < len(s):
            c = s[self.i]
            if c == "\\":
                if self.i + 1 < len(s) and s[self.i + 1] != "\n":
                    word = (word or "") + s[self.i + 1]
                self.i += 2
            elif c == "'":
                end = s.find("'", self.i + 1)
                if end < 0:
                    raise ValueError("閉じていない ' がある")
                word = (word or "") + s[self.i + 1 : end]
                self.i = end + 1
            elif c == '"':
                self.i += 1
                word = (word or "") + self._double_quoted()
            elif c in " \t\r":
                end_word()
                self.i += 1
            elif c == "\n":
                end_command()
                self.i += 1
                self._skip_heredoc_bodies()
            elif c == "#" and word is None:
                end = s.find("\n", self.i)
                self.i = len(s) if end < 0 else end
            elif c == "$" and s.startswith("$'", self.i):
                word = (word or "") + self._ansi_c_quoted()
            elif c == "$" and s.startswith("$((", self.i):
                self.i += 1
                word = (word or "") + self._arithmetic()
            elif c == "(" and s.startswith("((", self.i) and word is None and not words:
                word = self._arithmetic()
            elif c == "$" and s.startswith("$(", self.i):
                end_command()
                self.i += 2
                self.parse(close_paren=True)
            elif c == "(":
                end_command()
                depth += 1
                self.i += 1
            elif c == ")":
                end_command()
                self.i += 1
                if depth > 0:
                    depth -= 1
                elif close_paren:
                    return
            elif c == "`":
                end_command()
                self.i += 1
            elif c in "<>":
                end_word()
                if s.startswith("<<<", self.i):
                    self.i += 3
                    redirect_target = True
                elif s.startswith("<<", self.i):
                    self.i += 2
                    strip_tabs = s.startswith("-", self.i)
                    if strip_tabs:
                        self.i += 1
                    self.heredocs.append((self._delimiter(), strip_tabs))
                else:
                    self.i += 1
                    if self.i < len(s) and s[self.i] in "<>&|":
                        self.i += 1
                    redirect_target = True
            elif c in SEPARATORS:
                end_command()
                self.i += 1
            else:
                word = (word or "") + c
                self.i += 1
        if close_paren:
            raise ValueError("閉じていない $( がある")
        end_command()

    def _double_quoted(self):
        """`"` の次から閉じる `"` までを読み、中身を返す。中の `$(…)` はコマンドとして読む。"""
        s = self.s
        out = ""
        while self.i < len(s):
            c = s[self.i]
            if c == '"':
                self.i += 1
                return out
            if c == "\\" and self.i + 1 < len(s):
                if s[self.i + 1] in '$`"\\':
                    out += s[self.i + 1]
                elif s[self.i + 1] != "\n":
                    out += c + s[self.i + 1]
                self.i += 2
            elif s.startswith("$((", self.i):
                self.i += 1
                out += self._arithmetic()
            elif s.startswith("$(", self.i):
                self.i += 2
                self.parse(close_paren=True)
            else:
                out += c
                self.i += 1
        raise ValueError('閉じていない " がある')

    def _arithmetic(self):
        """`((` から対応する `))` までを読み、そのまま返す。

        中の `<<` はシフトで、ヒアドキュメントとして読むと後の行をすべて本文として飛ばす。中の置換は読まない。
        """
        s = self.s
        start = self.i
        depth = 0
        while self.i < len(s):
            if s[self.i] == "(":
                depth += 1
            elif s[self.i] == ")":
                depth -= 1
                if depth == 0:
                    self.i += 1
                    return s[start : self.i]
            self.i += 1
        raise ValueError("閉じていない (( がある")

    def _ansi_c_quoted(self):
        """`$'` から閉じる `'` までを読み、中身を返す（`\\` の後の 1 文字は引用を閉じない。エスケープは展開しない）。"""
        s = self.s
        start = self.i + 2
        self.i = start
        while self.i < len(s):
            if s[self.i] == "\\":
                self.i += 2
            elif s[self.i] == "'":
                self.i += 1
                return s[start : self.i - 1]
            else:
                self.i += 1
        raise ValueError("閉じていない $' がある")

    def _delimiter(self):
        """`<<` の後の区切りの語を、引用を外して読む。"""
        s = self.s
        while self.i < len(s) and s[self.i] in " \t":
            self.i += 1
        out = ""
        while self.i < len(s) and s[self.i] not in " \t\n;&|()<>":
            c = s[self.i]
            if c in "'\"":
                end = s.find(c, self.i + 1)
                if end < 0:
                    raise ValueError(f"閉じていない {c} がある")
                out += s[self.i + 1 : end]
                self.i = end + 1
            elif c == "\\" and self.i + 1 < len(s):
                out += s[self.i + 1]
                self.i += 2
            else:
                out += c
                self.i += 1
        if not out:
            raise ValueError("<< の後に区切りの語が無い")
        return out

    def _skip_heredoc_bodies(self):
        """改行の後で、待っているヒアドキュメントの本文を区切りの行まで飛ばす。"""
        s = self.s
        for delimiter, strip_tabs in self.heredocs:
            while self.i < len(s):
                end = s.find("\n", self.i)
                line = s[self.i : len(s) if end < 0 else end]
                self.i = len(s) if end < 0 else end + 1
                if (line.lstrip("\t") if strip_tabs else line) == delimiter:
                    break
        self.heredocs = []


def simple_commands(text):
    """コマンドの文字列を、実行される単純なコマンドの語の並びのリストにする。読めなければ ValueError。"""
    parser = _Parser(text)
    parser.parse()
    return parser.commands


def drop_options(argv, with_value=()):
    """先頭の `-` で始まる語（with_value の語は次の値も）を外す。`--` で終わる。

    短いオプションの束（`-iC`）は、束の最後の文字が値を取る短いオプションなら次の語も値として外す（束の途中なら、
    値は束の残りに付いている）。
    """
    short_with_value = {o[1] for o in with_value if len(o) == 2}
    while argv and argv[0].startswith("-") and argv[0] != "-":
        arg = argv[0]
        if arg == "--":
            return argv[1:]
        if arg in with_value:
            takes_next = True
        elif arg.startswith("--"):
            takes_next = False
        else:
            flags = arg[1:]
            k = next((k for k, f in enumerate(flags) if f in short_with_value), None)
            takes_next = k == len(flags) - 1
        argv = argv[2:] if takes_next else argv[1:]
    return argv


def _shell_command_string(args):
    """`bash` の引数から `-c` の文字列を返す。`-c` が無ければ空の並び。

    `-c` の文字列は、オプションの後の最初の語（bash(1) の `-c`）。`-o`、`-O`（`+o`、`+O`）と `--rcfile`、`--init-file` は
    次の語を値に取る。
    """
    command_string = False
    k = 0
    while k < len(args):
        arg = args[k]
        if arg == "--":
            k += 1
            break
        if arg in ("--rcfile", "--init-file"):
            k += 2
        elif arg.startswith("--"):
            k += 1
        elif arg[:1] in "-+" and len(arg) > 1:
            command_string = command_string or "c" in arg[1:]
            k += 2 if arg[-1] in "oO" else 1
        else:
            break
    return args[k] if command_string and k < len(args) else []


def unwrap(argv):
    """先頭の代入と包むコマンドを外した語の並びを返す。中のコマンドが文字列（bash -c）なら、その文字列を返す。

    包むコマンドの下に実行するコマンドが無ければ（`nix develop` だけ、`bash` だけ）、空の並びを返す。
    """
    while argv:
        if ASSIGNMENT.match(argv[0]) or argv[0] in RESERVED:
            argv = argv[1:]
            continue
        name = os.path.basename(argv[0])
        rest = argv[1:]
        if name in ("command", "exec", "builtin", "nohup", "time", "sudo"):
            argv = drop_options(rest)
        elif name == "env":
            argv = drop_options(rest, with_value=("-u", "-C", "--unset", "--chdir"))
        elif name == "nice":
            argv = drop_options(rest, with_value=("-n",))
        elif name == "timeout":
            argv = drop_options(
                rest, with_value=("-s", "-k", "--signal", "--kill-after")
            )[1:]
        elif name == "taskset":
            argv = drop_options(rest)[1:]
        elif name == "nix" and rest[:1] in (["develop"], ["shell"]):
            k = next((k for k, a in enumerate(rest) if a in ("-c", "--command")), None)
            if k is None:
                return []
            argv = drop_options(rest[k + 1 :])
        elif name in SHELLS:
            return _shell_command_string(rest)
        elif PYTHONS.fullmatch(name) and rest[:2] == ["-m", "dvc"]:
            return ["dvc", *rest[2:]]
        else:
            return argv
    return argv


def _runs_dvc_command(argv, subcommand):
    unwrapped = unwrap(argv)
    if isinstance(unwrapped, str):  # bash -c '<文字列>'
        return runs_dvc(unwrapped, subcommand)
    if not unwrapped or os.path.basename(unwrapped[0]) != "dvc":
        return False
    sub = drop_options(unwrapped[1:], with_value=("--cd",))
    return sub[:1] == [subcommand]


def runs_dvc(command, subcommand):
    """コマンドが dvc の subcommand（checkout、repro など）を直接呼ぶか。コマンドを読めなければ ValueError。"""
    return any(_runs_dvc_command(argv, subcommand) for argv in simple_commands(command))
