#!/usr/bin/env python3
"""
Print a Mermaid diagram of which analysis modules import which, from the actual import
statements (only the analysis' own .py files; pocket-coffea, numpy, ... are left out).

    python3 ../common/module_graph.py . ../common > module_graph.mmd
    python3 ../common/module_graph.py . ../common --from run_pipeline makeDatacard

Paste the output into a Markdown file inside ```mermaid fences (Obsidian, GitHub and VS Code
render it), or into https://mermaid.live. With --from, only the modules reachable from the given
entry points are drawn.
"""
import argparse
import ast
import os


def local_modules(dirs):
    """{module name: path} of the .py files in the given directories."""
    modules = {}
    for d in dirs:
        for name in sorted(os.listdir(d)):
            if name.endswith(".py") and not name.startswith("_"):
                modules.setdefault(name[:-3], os.path.join(d, name))
    return modules


def imports_of(path, modules):
    tree = ast.parse(open(path).read(), filename=path)
    found = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            found.update(a.name.split(".")[0] for a in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module and node.level == 0:
            found.add(node.module.split(".")[0])
    return sorted(m for m in found if m in modules)


def reachable(graph, starts):
    seen, todo = set(), list(starts)
    while todo:
        m = todo.pop()
        if m not in seen:
            seen.add(m)
            todo += graph.get(m, [])
    return seen


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dirs", nargs="+", help="directories with the analysis modules")
    ap.add_argument("--from", dest="starts", nargs="*", help="only modules reachable from these")
    args = ap.parse_args()

    modules = local_modules(args.dirs)
    graph = {m: [i for i in imports_of(p, modules) if i != m] for m, p in modules.items()}
    keep = reachable(graph, args.starts) if args.starts else set(graph)

    print("flowchart LR")
    for m in sorted(keep):
        for dep in graph[m]:
            if dep in keep:
                print(f"    {m} --> {dep}")
    for m in sorted(keep):
        if not graph[m] and not any(m in graph[o] for o in keep):
            print(f"    {m}")


if __name__ == "__main__":
    main()
