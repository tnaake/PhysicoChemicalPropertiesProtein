## sample sequences taken from
## https://www.bioinformatics.org/sms2/protein_gravy.html
test_that("checkAA works.", {
    
    expect_error(checkAA(1),
                 "'aa' has to be an 'AAString' object")
    expect_error(checkAA(c("MQ", "AELS")),
                 "'aa' has to be an 'AAString' object")
    seq_aa <- AAString("xycAME")
    expect_error(checkAA(seq_aa), 
        "'aa' contains letters not defined in the amino acid alphabet")
    seq_aa <- AAString("MQKSPLEKASFISKLFFSWTTPILRKGYR")
    expect_equal(checkAA(seq_aa), "MQKSPLEKASFISKLFFSWTTPILRKGYR")
})
