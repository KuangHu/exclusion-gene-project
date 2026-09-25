Subject: lr6-OST006a down since 18 Sep — and the default layout is still allocating onto it

Hello,

I'm reporting a Lustre OST failure on /global/scratch that has been ongoing since
18 September and is still present today. There are two separate impacts, and I
think the second one may not be on your radar.

## 1. Reads from lr6-OST006a fail deterministically

Any file whose first stripe component landed on lr6-OST006a (index 106) is
unreadable. This is not intermittent: the same files fail on every attempt, and
files on other OSTs are unaffected.

`lfs check osts` does report the fault:

```
$ lfs check osts
lfs check: error: check 'lr6-OST006a-osc-ffff97dd43959000':
           Cannot send after transport endpoint shutdown (108)
```

`lfs df`, however, simply omits the OST — it lists 29 targets where
`lfs check osts` lists 30 — so a capacity or count-based health check looks
normal unless you already know the expected total.

Reproducer:

```
$ head -c 32 /global/scratch/users/kh36969/exclusion_gene/cds_cache_full/cds_002.faa
head: error reading '.../cds_002.faa': Input/output error

$ lfs getstripe /global/scratch/users/kh36969/exclusion_gene/cds_cache_full/cds_002.faa
  lcm_layout_gen:    4
  lcm_mirror_count:  1
  lcm_entry_count:   3
    lcme_id:             1
    lcme_extent.e_start: 0
    lcme_extent.e_end:   65536
      lmm_stripe_count:  1
      lmm_stripe_size:   65536
      lmm_stripe_offset: 106
      lmm_pool:          ddn_nvme2
      lmm_objects:
      - 0: { l_ost_idx: 106, l_fid: [0x6400069a4:0x1467ae6:0x0] }
    [components 2 and 3 are on healthy OSTs]
```

The error surfaces as EIO on read and as ESHUTDOWN on stat, so `ls -l` on an
affected directory prints `-????????? ? ? ?` rows alongside normal ones.

Affected files under my account (kh36969):

```
$ lfs find /global/scratch/users/kh36969 --ost 106 | wc -l
1982920
```

Just under 2 million files, which is consistent with the ~25% allocation rate in
point 2 below. That count is files holding an object on OST 106 anywhere in their
layout; a file is wholly unreadable when the object is in component 1, and
partially unreadable when it is in a later component. By directory, for the
datasets I actively depend on:

```
funcannot_dbs   21,679      reference databases, HMM profile sets
exclusion_gene     107      derived caches and models
ice_db               7
plsdb_dedup          7
plsdb                3      including sequences.fasta, nuccore.csv, taxonomy.csv
```

## 2. New writes are still being allocated onto the failed OST

The default PFL layout places a new file's first component in pool `ddn_nvme2`,
and lr6-OST006a is one of that pool's four members. New files are therefore
still being created on the dead OST, and they are unreadable from the moment
they are written — with no error raised at write time.

Measured directly, writing 24 fresh 1 MB files into a default-layout directory:

```
first component allocated on OST 106:   6 / 24
readable immediately after writing:    18 / 24
```

The six unreadable files are exactly the six that landed on 106.

I hit this in production before I measured it: a decompression job I launched
wrote its output straight onto OST 106. I cancelled the job once I checked the
output file's layout, but had I not checked, it would have run to completion and
reported success while producing a file nothing could read.

If this is correct, the impact is not limited to users reading older data. Any
user writing with the default layout is currently producing unreadable output
silently, at roughly a 25% rate.

## Questions

1. Is the data on lr6-OST006a expected to be recoverable once the OST is brought
   back, or should we treat it as lost? I have regeneration paths for most of my
   affected files but would rather not spend the compute if the originals are
   coming back.
2. Is there an ETA for the OST, or for a decision on its status?
3. In the meantime, is pinning new writes away from `ddn_nvme2` the recommended
   workaround, or is there something better? I have pinned my output directories
   with `lfs setstripe -p lr6.ddn_hdd` (22 OSTs, none of them 106) and verified
   that probe files read back correctly.

Given point 2, a broader notice to users on lr6 may be warranted.

Happy to provide additional affected paths, `lfs getstripe` output, or job IDs.

Thanks,
Kuangh
kh36969
