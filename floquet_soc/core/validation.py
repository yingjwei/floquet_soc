# -*- coding: utf-8 -*-
"""Minimal validation framework: PASS / FAIL / WARN / BLOCKED / SKIP reports."""
from __future__ import annotations

import json
import os
from dataclasses import dataclass, field, asdict

PASS, FAIL, WARN, BLOCKED, SKIP = "PASS", "FAIL", "WARN", "BLOCKED", "SKIP"


@dataclass
class Check:
    name: str
    status: str
    detail: str = ""
    value: float | None = None
    tol: float | None = None

    def __post_init__(self):
        if self.status not in (PASS, FAIL, WARN, BLOCKED, SKIP):
            raise ValueError(f"bad status {self.status!r}")


@dataclass
class ValidationReport:
    title: str
    checks: list = field(default_factory=list)

    def add(self, name, status, detail="", value=None, tol=None):
        self.checks.append(Check(name, status, detail, value, tol))
        return self

    def extend(self, items):
        for it in items:
            if isinstance(it, Check):
                self.checks.append(it)
            else:
                self.add(*it)
        return self

    def count(self, status):
        return sum(1 for c in self.checks if c.status == status)

    @property
    def ok(self):
        return self.count(FAIL) == 0 and self.count(BLOCKED) == 0

    def summary(self):
        return {s: self.count(s) for s in (PASS, FAIL, WARN, BLOCKED, SKIP)}

    def to_markdown(self) -> str:
        lines = [f"# {self.title}", "",
                 f"Summary: " + ", ".join(
                     f"**{k}**={v}" for k, v in self.summary().items()), "",
                 "| check | status | detail |", "|---|---|---|"]
        for c in self.checks:
            lines.append(f"| `{c.name}` | **{c.status}** | {c.detail} |")
        return "\n".join(lines) + "\n"

    def to_json(self) -> str:
        return json.dumps(dict(title=self.title, summary=self.summary(),
                               checks=[asdict(c) for c in self.checks]),
                          ensure_ascii=False, indent=2)

    def save(self, outdir: str, stem: str):
        os.makedirs(outdir, exist_ok=True)
        md = os.path.join(outdir, stem + ".md")
        js = os.path.join(outdir, stem + ".json")
        with open(md, "w", encoding="utf-8") as f:
            f.write(self.to_markdown())
        with open(js, "w", encoding="utf-8") as f:
            f.write(self.to_json())
        return md, js

    def print(self):
        print(f"\n== {self.title} ==")
        for c in self.checks:
            mark = {"PASS": "ok  ", "FAIL": "FAIL", "WARN": "warn",
                    "BLOCKED": "BLOCK", "SKIP": "skip"}[c.status]
            print(f"  [{mark}] {c.name}: {c.detail}")
        print("  summary:", self.summary())


def check_bool(name, cond, detail, value=None, tol=None):
    return Check(name, PASS if cond else FAIL, detail, value, tol)
