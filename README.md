# statest

Some statistical tests

# Data oil_uv
With the test design `~uv + oil + uv:oil`, we can test the effect of uv on oil, and also taking the different oil does into consideration. I think both DESeq2 and edgeR should be able to do the test. But since I found the exactly same example from edgeR’s tutorial ([edgeR](https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf) Section 3.3.4), I conducted the test using edgeR.

I’ve done a factorial test to see the effect of UV light on oil (section `effect of uv on oil` in the script and output file `dea_uv_effect_on_oil.tsv`). There are two columns of fold change: `logFC.uvyes.oilhigh` and `logFC.uvyes.oillow`, meaning the following contrast:

`(uvyes.oilhigh - uvyes.oilctrl) - (uvno.oilhigh - uvno.oilctrl)`

`(uvyes.oillow - uvyes.oilctrl) - (uvno.oillow - uvno.oilctrl)`

Considering “It is known that UV light alters the chemical structure of PAHs and the “modified” PAHs compounds can be more toxic”, I’ve done another two tests: oil’s effect when there is uv light, and oil’s effect when there is not uv light. There are two output files `dea_oil_effect_with_uv.tsv` and `dea_oil_effect_without_uv.tsv`. There are two columns of fold change in each of them: `logFC.oilhigh` and `logFC.oillow`. 

These two FCs in `with_uv` means the following contrast:

`uvyes.oilhigh - uvyes.oilctrl`

`uvyes.oillow - uvyes.oilctrl`

These two FCs in `without_uv` means the following contrast:

`uvno.oilhigh - uvno.oilctrl`

`uvno.oillow - uvno.oilctrl`

Comparing these two output files, it is confirmed that ““modified” PAHs compounds can be more toxic”: more differentially expressed genes and higher FC (but of course not for every gene) in the `with_uv` condition.

Update on July 6, 2020:
script: script/simple_test_with_dose.R
Doing simple tests:
- with uv: low vs. ctrl, high vs. ctrl
- without uv: low vs. ctrl, high vs. ctrl



Update on Sept. 22th, 2021:

Script: test_uv.R

To test the effect of UV without oil: 

- UV_yes_oil_ctrl vs. UV_no_oil_ctrl