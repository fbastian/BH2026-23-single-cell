#!/usr/bin/env python
"""Font-size audit of figure PDFs against the Nature specifications: all text 5-7 pt, panel letters 8 pt bold
(single lowercase letters at 8 pt are accepted as panel letters). Reports every text span outside 5-8 pt, per file.
Page size (corrected 2026-11-06, user-pasted Nature spec): width <= 183 mm and height <= 170 mm (247 mm = full page depth,
reference only; the old <= 247 mm check was wrong). --working: pages that are review copies, size reported as WORKING not FAIL.
Only meaningful for PDFs saved at their print size (all ac_* Nature figures are).

    masld_atlas env:  python -I code/annotation_check/check_nature_fonts.py <pdf> [<pdf> ...]
"""
import sys
from collections import Counter

import pymupdf


def audit(path):
    sizes, bad = Counter(), []
    doc = pymupdf.open(path)
    for page in doc:
        for b in page.get_text('dict')['blocks']:
            for line in b.get('lines', []):
                for sp in line['spans']:
                    t, s = sp['text'].strip(), round(sp['size'], 1)
                    if not t:
                        continue
                    sizes[s] += 1
                    if s < 4.95 or s > 8.05 or (s > 7.05 and not (len(t) == 1 and t.islower())):
                        bad.append((s, t[:40]))
    w, h = doc[0].rect.width / 72 * 25.4, doc[0].rect.height / 72 * 25.4
    return sizes, bad, (w, h)


WORKING = '--working' in sys.argv[1:]


def main():
    fail = 0
    for p in [a for a in sys.argv[1:] if a != '--working']:
        sizes, bad, (w, h) = audit(p)
        size_ok = w <= 183.5 and h <= 170.5
        ok = not bad and (size_ok or WORKING)
        fail += not ok
        tag = 'OK  ' if (ok and size_ok) else ('WORKING (size not Nature)' if ok else 'FAIL')
        print(f'{tag} {p}  {w:.0f} x {h:.0f} mm  sizes {dict(sorted(sizes.items()))}' + ('' if size_ok else '  [page > 183 x 170 mm]'))
        for s, t in sorted(set(bad))[:12]:
            print(f'      {s} pt: {t!r}')
    sys.exit(1 if fail else 0)


if __name__ == '__main__':
    main()
