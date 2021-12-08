test_that("calculateIsoelectricPoint", {
    
    expect_error(calculateIsoelectricPoint(seq, "foo"), 
                 "'arg' should be one of")
    
    ## values taken from http://isoelectric.org/calculate.php
    ## P99029-1, experimental isoelectric points (different sources): 
    ## 7.84, 7.65, 7.54
    seq <- "MLQLGLRVLGCKASSVLRASTCLAGRAGRKEAGWECGGARSFSSSAVTMAPIKVGDAIPSVEVFEGEPGKKVNLAELFKGKKGVLFGVPGAFTPGCSKTHLPGFVEQAGALKAKGAQVVACLSVNDVFVIEEWGRAHQAEGKVRLLADPTGAFGKATDLL
LDDSLVSLFGNRRLKRFSMVIDNGIVKALNVEPDGTGLTCSLAPNILSQL"
    expect_equal(calculateIsoelectricPoint(seq, "EMBOSS"), 9.147, 
                 tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "DTASelect"), 8.848, 
                 tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Solomon"), 9.101,
                 tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Sillero"), 9.291, 
                 tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Rodwell"), 9.03,
                 tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Lehninger"), 9.127,
                 tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Toseland"), 8.18, 
                 tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Thurlkill"), 9.02, 
                 tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Nozaki"), 9.542, 
                 tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "IPC_protein"), 8.054,
                 tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "IPC_peptide"), 9.101, 
                 tolerance = 5e-02)
    
    ## P40185-1 experimental isoelectric point 6.89
    seq <- "MFLRNSVLRTAPVLRRGITTLTPVSTKLAPPAAASYSQAMKANNFVYVSGQIPYTPDNKPVQGSISEKAEQVFQNVKNILAESNSSLDNIVKVNVFLADMKNFAEFNSVYAKHFHTHKPARSCVGVASLPLNVDLEMEVIAVEKN"
    expect_equal(calculateIsoelectricPoint(seq, "EMBOSS"), 9.764, 
                 tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "DTASelect"), 9.269, 
                 tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Solomon"), 9.694, 
                 tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Sillero"), 9.542, 
                 tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Rodwell"), 9.918, 
                 tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Lehninger"), 9.667, 
                 tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Toseland"), 9.357, 
                 tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Thurlkill"), 9.443, 
                 tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Nozaki"), 9.421, 
                 tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "IPC_protein"), 8.64, 
                 tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "IPC_peptide"), 9.687, 
                 tolerance = 5e-02)
    
    ## Q9XFT3-1 experimental isoelectric point 8.07
    seq <- "MASMGGLHGASPAVLEGSLKINGSSRLNGSGRVAVAQRSRLVVRAQQSEETSRRSVIGLVAAGLAGGSFVQAVLADAISIKVGPPPAPSGGLPAGTDNSDQARDFALALKDRFYLQPLPPTEAAARAKESAKDIINVKPLIDRKAWPYVQNDLRSKASYLRYDLNTIISSKPKDEKKSLKDLTTKLFDTIDNLDYAAKKKSPSQAEKYYAETVSALNEVLAKLG"
    expect_equal(calculateIsoelectricPoint(seq, "EMBOSS"), 10.222, 
                 tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "DTASelect"), 9.645, 
                 tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Solomon"), 10.046, 
                 tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Sillero"), 9.921, 
                 tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Rodwell"), 10.498, 
                 tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Lehninger"), 10.017, 
                 tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Toseland"), 9.824,
                 tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Thurlkill"), 9.866, 
                 tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Nozaki"), 9.784,
                 tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "IPC_protein"), 8.958, 
                 tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "IPC_peptide"), 10.046, 
                 tolerance = 5e-02)
    
    ## Q9BVP2-1 experimental isoelectric point 6.47
    seq <- "MKRPKLKKASKRMTCHKRYKIQKKVREHHRKLRKEAKKRGHKKPRKDPGVPNSAPFKEALLREAELRKQRLEELKQQQKLDRQKELEKKRKLETNPDIKPSNVEPMEKEFGLCKTENKAKSGKQNSKKLYCQELKKVIEASDVVLEVLDARDPLGCRCPQVEEAIVQSGQKKLVLILNKSDLVPKENLESWLNYLKKELPTVVFRASTKPKDKGKITKRVKAKKNAAPFRSEVCFGKEGLWKLLGGFQETCSKAIRVGVIGFPNVGKSSIINSLKQEQMCNVGVSMGLTRSMQVVPLDKQITIIDSPSFIVSPLNSSSALALRSPASIEVVKPMEAASAILSQADARQVVLKYTVPGYRNSLEFFTVLAQRRGMHQKGGIPNVEGAAKLLWSEWTGASLAYYCHPPTSWTPPPYFNESIVVDMKSGFNLEELEKNNAQSIRAIKGPHLANSILFQSSGLTNGIIEEKDIHEELPKRKERKQEEREDDKDSDQETVDEEVDENSSGMFAAEETGEALSEETTAGEQSTRSFILDKIIEEDDAYDFSTDYV"
    expect_equal(calculateIsoelectricPoint(seq, "EMBOSS"), 9.41,
                 tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "DTASelect"), 8.933, 
                 tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Solomon"), 9.249, 
                 tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Sillero"), 9.304,
                 tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Rodwell"), 9.579, 
                 tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Lehninger"), 9.237, 
                 tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Toseland"), 8.903, 
                 tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Thurlkill"), 9.149,
                 tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Nozaki"), 9.363, 
                 tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "IPC_protein"), 8.081, 
                 tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "IPC_peptide"), 9.251, 
                 tolerance = 5e-02)
})
