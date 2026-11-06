"""Check the agent guidance files: relative links, heading anchors and named tests.

    python tools/check_guidance.py

Scans AGENTS.md, ai_skills/, docs/lessons/, docs/pivots.md and examples/README.md. For each relative Markdown link it
checks that the target file exists and, when the link has a #fragment, that the target page has a heading or an explicit
`<a id="...">` with that name. Heading anchors follow GitHub's slugger: lower case, punctuation removed (hyphens and
underscores kept), each space turned into a hyphen, and -1, -2 ... appended to repeated anchors. Text inside fenced code
blocks is ignored. For each `tests/<file>.py::<name>` reference it checks that the test function exists. External links
are not fetched. Exit status 1 on any problem; the test suite runs it (tests/test_guidance_links.py).
"""
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCOPE = ['AGENTS.md', 'ai_skills', 'docs/lessons', 'docs/pivots.md', 'examples/README.md']

LINK = re.compile(r'\[[^\]]*\]\(([^)\s]+)\)')
TESTREF = re.compile(r'(tests/[\w/]+\.py)::(\w+)')
HEADING = re.compile(r'^ {0,3}#{1,6}[ \t]+(.*?)(?:[ \t]+#+)?[ \t]*$')
EXPLICIT_ID = re.compile(r'<a\s+(?:id|name)="([^"]+)"')
FENCE = re.compile(r'^ {0,3}(`{3,}|~{3,})')


def slug(heading):
    """GitHub's anchor for one heading's text (before de-duplication)."""
    text = re.sub(r'!?\[([^\]]*)\]\([^)]*\)', r'\1', heading)   # a link renders as its text
    text = text.lower()
    text = re.sub(r'[^\w\- ]', '', text)                         # drop punctuation; keep - and _
    return text.replace(' ', '-')


def prose_lines(text):
    """Yield the lines outside fenced code blocks (``` or ~~~, closed by the same character, at least as long)."""
    fence = None
    for line in text.splitlines():
        m = FENCE.match(line)
        if fence is None:
            if m:
                fence = m.group(1)
                continue
            yield line
        elif m and m.group(1)[0] == fence[0] and len(m.group(1)) >= len(fence) and not line.strip()[len(m.group(1)):]:
            fence = None


def anchors(path, cache={}):
    if path not in cache:
        found, counts = set(), {}
        for line in prose_lines(open(path, encoding='utf-8').read()):
            found.update(EXPLICIT_ID.findall(line))
            m = HEADING.match(line)
            if m:
                base = result = slug(m.group(1))
                while result in counts:
                    counts[base] += 1
                    result = '%s-%d' % (base, counts[base])
                counts[result] = 0
                found.add(result)
        cache[path] = found
    return cache[path]


def files(root):
    for entry in SCOPE:
        p = os.path.join(root, entry)
        if os.path.isfile(p):
            yield p
        elif os.path.isdir(p):
            for d, _, names in os.walk(p):
                for n in sorted(names):
                    if n.endswith('.md'):
                        yield os.path.join(d, n)


def check(root=ROOT):
    """Return (number of links checked, list of problems) for the guidance files under root."""
    problems, nlinks = [], 0
    for src in files(root):
        rel = os.path.relpath(src, root)
        text = open(src, encoding='utf-8').read()
        prose = '\n'.join(prose_lines(text))
        for target in LINK.findall(prose):
            if re.match(r'[a-z]+:', target):
                continue
            nlinks += 1
            path, _, frag = target.partition('#')
            dest = os.path.normpath(os.path.join(os.path.dirname(src), path)) if path else src
            if not os.path.exists(dest):
                problems.append('%s: missing file %s' % (rel, target))
            elif frag and dest.endswith('.md') and frag not in anchors(dest):
                problems.append('%s: missing anchor %s' % (rel, target))
        for tfile, name in TESTREF.findall(text):
            tpath = os.path.join(root, tfile)
            if not os.path.exists(tpath):
                problems.append('%s: missing test file %s' % (rel, tfile))
            elif not re.search(r'^\s*def %s\(' % re.escape(name), open(tpath, encoding='utf-8').read(), flags=re.M):
                problems.append('%s: missing test %s::%s' % (rel, tfile, name))
    return nlinks, problems


def main():
    nlinks, problems = check()
    for p in problems:
        print(p)
    print('%d links checked, %d problems' % (nlinks, len(problems)))
    return 1 if problems else 0


if __name__ == '__main__':
    sys.exit(main())
