# hapfuse

This is my own fork of the SNPTools::hapfuse tool in [SNPTools v1.0](http://sourceforge.net/projects/snptools).

## Differences to SNPTools::hapfuse

* **Dependency on htslib for parsing input chunks**.
  This allows me to avoid nasty buffer overruns when input VCFs contain too many samples. Also allows for the use of BCF files.

* **Input VCFs must contain the GT, and APP or GP field**.
  APP stands for phred-scaled allelic probability.  SNPTools uses the AP field, which stands for allelic probability.

* **Input VCFs can have fields with arbitrary precision**.  SNPTools::hapfuse uses the AP field with fixed precision floats while hapfuse v0.4 accepts APP and GP fields of arbitrary precision.  

* **Speed and Reliability**.
  Hapfuse is now way faster than SNPTools::hapfuse as it uses htslib for VCF/BCF parsing.  Hapfuse also concurrently loads, fuses and writes chunks, so expect to see CPU usage at about %200.

* **Regression Tests**. Hapfuse comes with regression tests to make sure it works!

* **Concurrency**. Hapfuse now uses concurrency to simultaneously load, merge and write VCF/BCF files.  This gives about a 2 fold speedup.

## WARNING
Hapfuse does not work with data generated by SNPTools::impute.  This is because I do not use the AP field and am too lazy to implement support for it.  Patches welcome.