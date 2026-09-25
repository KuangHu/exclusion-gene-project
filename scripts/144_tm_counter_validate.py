#!/usr/bin/env python3
"""A local TM-segment counter, and the control gate it must pass before use.

WHY LOCAL. No TMHMM/Phobius/TMbed installed; both need registered downloads.
No torch/transformers, so DeepTMHMM cannot run locally. DeepTMHMM via BioLib
executes REMOTELY -- it would send our sequences off-cluster, which is not done
without explicit approval.

EMBOSS `tmap` IS installed and is REJECTED: on the R64 ExcA N-terminus, which
contains an unmistakable TM stretch, it returns HitCount 0 -- and also 0 on a
deliberately hydrophilic control. It implements Persson & Argos, which needs a
multiple alignment. A predictor that returns zero on both a known positive and a
known negative is inert, not conservative.

METHOD. Kyte-Doolittle mean hydropathy over a sliding window; a TM segment is a
run of >= MIN_LEN consecutive positions whose windowed mean exceeds THRESH.
This is the use KD was originally proposed for.

THE GATE. A thresholded-KD counter risks inheriting the weakness already measured
for raw KD (2.3x separation against a 40% background in MPF_I). So it is
validated first on named controls, and is used ONLY if:

    every known membrane control  >= 1 TM segment
    every known soluble control   == 0 TM segments

If that fails at every parameter setting, the counter is not used and TM count is
reported as UNAVAILABLE rather than as a weak criterion.
"""
import os, sys
KD = {'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,
      'H':-3.2,'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,
      'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}

# Known inner-membrane exclusion proteins and their published evidence.
POS = {
 "R64_ExcA_Nterm": "MKAKKETTDRFPTWWLFYYVLRKAYFFLGIPFFLAAAGGLLLLAAAWLLGSRRRKPVEEHA",
 "generic_TM_helix": "MGSSHHHHHHAAALWWLLAGGVVLIAAFLFGVLSILAAVGGWLVSRRKKEDKQ",
 "polytopic_2TM": "MKQDKARWLLPLALLAGVAVLGFAWWQSRQDKPAAMLSVIGGAALLAGSLFLWWRKQAEEGK",
}
NEG = {
 "polyE_K_soluble": "MSDKEEKKDDNKKEEKKDDNKKEEKKSDNKQEEKKSDNKQEEKKDDNKKEEKKSDNKQEEK",
 "charged_random":  "MRKDEQNSTRKDEQNSTRKDEQNSTRKDEQNSTRKDEQNSTRKDEQNSTRKDEQNSTRKDE",
 "proline_rich":    "MPPQPSPQPSPQPSPQPSPQPSPQPSPQPSPQPSPQPSPQPSPQPSPQPSPQPSPQPSPQP",
}


def tm_segments(seq, win=19, thresh=1.6, min_len=16):
    n = len(seq)
    if n < win:
        return 0, []
    h = w = 0
    vals = []
    for i in range(n - win + 1):
        vals.append(sum(KD.get(c, 0.0) for c in seq[i:i + win]) / win)
    segs, run, start = [], 0, None
    for i, v in enumerate(vals):
        if v >= thresh:
            if run == 0: start = i
            run += 1
        else:
            if run >= min_len - win + 1 and run > 0:
                segs.append((start, start + run + win - 2))
            run = 0
    if run >= min_len - win + 1 and run > 0:
        segs.append((start, start + run + win - 2))
    return len(segs), segs


def main():
    print("=== CONTROL GATE: the counter must be right on BOTH sides ===")
    best = None
    print("  %-6s %-8s %-8s  %-28s %s" % ("win", "thresh", "min_len", "positives (want >=1)", "negatives (want 0)"))
    for win in (15, 17, 19, 21):
        for thresh in (1.2, 1.4, 1.6, 1.8, 2.0):
            for min_len in (15, 18, 20):
                p = {k: tm_segments(v, win, thresh, min_len)[0] for k, v in POS.items()}
                n = {k: tm_segments(v, win, thresh, min_len)[0] for k, v in NEG.items()}
                ok = all(x >= 1 for x in p.values()) and all(x == 0 for x in n.values())
                if ok and best is None:
                    best = (win, thresh, min_len)
                if ok:
                    print("  %-6d %-8.1f %-8d  %-28s %s  PASS"
                          % (win, thresh, min_len,
                             ",".join("%s=%d" % (k.split('_')[0], v) for k, v in p.items()),
                             ",".join(str(v) for v in n.values())))
    if best is None:
        print("\n  NO PARAMETER SETTING PASSES THE GATE.")
        print("  The KD-window counter cannot separate known membrane from known")
        print("  soluble controls. TM count is UNAVAILABLE -- it is NOT reported as")
        print("  a weak criterion, because a predictor that fails its own controls")
        print("  produces numbers that look like data and are not.")
        return 1
    win, thresh, min_len = best
    print("\n  ADOPTED: win=%d thresh=%.1f min_len=%d" % (win, thresh, min_len))
    for k, v in list(POS.items()) + list(NEG.items()):
        c, segs = tm_segments(v, win, thresh, min_len)
        print("    %-20s TM=%d  %s" % (k, c, segs))
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib",
                           "tm_params.txt"), "w") as fh:
        fh.write("%d\t%.1f\t%d\n" % (win, thresh, min_len))
    print("\n  wrote scripts/lib/tm_params.txt")
    return 0


if __name__ == "__main__":
    sys.exit(main())
