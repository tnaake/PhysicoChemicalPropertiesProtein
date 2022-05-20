test_that("calculateIsoelectricPoint works.", {
    
    ## check method
    seq_aa <- AAString("TEST")
    expect_error(calculateIsoelectricPoint(seq_aa, method = "foo"), 
                 "'arg' should be one of")
    
    ## check letters in amino acids
    seq_aa <- AAString("aaa")
    expect_error(calculateIsoelectricPoint(seq_aa, method = "EMBOSS"), 
                 "'aa' contains letters not defined in the amino acid alphabet")
    
    ## values taken from http://isoelectric.org/calculate.php
    ## P99029-1, experimental isoelectric points (different sources): 
    ## 7.84, 7.65, 7.54
    seq_aa <- AAString("MLQLGLRVLGCKASSVLRASTCLAGRAGRKEAGWECGGARSFSSSAVTMAPIKVGDAIPSVEVFEGEPGKKVNLAELFKGKKGVLFGVPGAFTPGCSKTHLPGFVEQAGALKAKGAQVVACLSVNDVFVIEEWGRAHQAEGKVRLLADPTGAFGKATDLLLDDSLVSLFGNRRLKRFSMVIDNGIVKALNVEPDGTGLTCSLAPNILSQL")
    expect_equal(calculateIsoelectricPoint(seq_aa, "EMBOSS"), 9.147, 
                 tolerance = 1e-03)
    expect_equal(calculateIsoelectricPoint(seq_aa, "DTASelect"), 8.848, 
                 tolerance = 1e-03)
    expect_equal(calculateIsoelectricPoint(seq_aa, "Solomon"), 9.101,
                 tolerance = 1e-03)
    expect_equal(calculateIsoelectricPoint(seq_aa, "Sillero"), 9.291, 
                 tolerance = 1e-03)
    expect_equal(calculateIsoelectricPoint(seq_aa, "Rodwell"), 9.03,
                 tolerance = 1e-03)
    expect_equal(calculateIsoelectricPoint(seq_aa, "Lehninger"), 9.127,
                 tolerance = 1e-03)
    expect_equal(calculateIsoelectricPoint(seq_aa, "Toseland"), 8.18, 
                 tolerance = 1e-03)
    expect_equal(calculateIsoelectricPoint(seq_aa, "Thurlkill"), 9.02, 
                 tolerance = 1e-03)
    expect_equal(calculateIsoelectricPoint(seq_aa, "Nozaki"), 9.542, 
                 tolerance = 1e-03)
    expect_equal(calculateIsoelectricPoint(seq_aa, "IPC_protein"), 8.054,
                 tolerance = 1e-03)
    expect_equal(calculateIsoelectricPoint(seq_aa, "IPC_peptide"), 9.101, 
                 tolerance = 1e-03)
    
    ## P40185-1 experimental isoelectric point 6.89
    seq_aa <- AAString("MFLRNSVLRTAPVLRRGITTLTPVSTKLAPPAAASYSQAMKANNFVYVSGQIPYTPDNKPVQGSISEKAEQVFQNVKNILAESNSSLDNIVKVNVFLADMKNFAEFNSVYAKHFHTHKPARSCVGVASLPLNVDLEMEVIAVEKN")
    expect_equal(calculateIsoelectricPoint(seq_aa, "EMBOSS"), 9.764, 
                 tolerance = 1e-03)
    expect_equal(calculateIsoelectricPoint(seq_aa, "DTASelect"), 9.269, 
                 tolerance = 1e-03)
    expect_equal(calculateIsoelectricPoint(seq_aa, "Solomon"), 9.694, 
                 tolerance = 1e-03)
    expect_equal(calculateIsoelectricPoint(seq_aa, "Sillero"), 9.542, 
                 tolerance = 1e-03)
    expect_equal(calculateIsoelectricPoint(seq_aa, "Rodwell"), 9.918, 
                 tolerance = 1e-03)
    expect_equal(calculateIsoelectricPoint(seq_aa, "Lehninger"), 9.667, 
                 tolerance = 1e-03)
    expect_equal(calculateIsoelectricPoint(seq_aa, "Toseland"), 9.357, 
                 tolerance = 1e-03)
    expect_equal(calculateIsoelectricPoint(seq_aa, "Thurlkill"), 9.443, 
                 tolerance = 1e-03)
    expect_equal(calculateIsoelectricPoint(seq_aa, "Nozaki"), 9.421, 
                 tolerance = 1e-03)
    expect_equal(calculateIsoelectricPoint(seq_aa, "IPC_protein"), 8.64, 
                 tolerance = 1e-03)
    expect_equal(calculateIsoelectricPoint(seq_aa, "IPC_peptide"), 9.687, 
                 tolerance = 1e-03)
    
    ## Q9XFT3-1 experimental isoelectric point 8.07
    seq_aa <- AAString("MASMGGLHGASPAVLEGSLKINGSSRLNGSGRVAVAQRSRLVVRAQQSEETSRRSVIGLVAAGLAGGSFVQAVLADAISIKVGPPPAPSGGLPAGTDNSDQARDFALALKDRFYLQPLPPTEAAARAKESAKDIINVKPLIDRKAWPYVQNDLRSKASYLRYDLNTIISSKPKDEKKSLKDLTTKLFDTIDNLDYAAKKKSPSQAEKYYAETVSALNEVLAKLG")
    expect_equal(calculateIsoelectricPoint(seq_aa, "EMBOSS"), 10.222, 
                 tolerance = 1e-03)
    expect_equal(calculateIsoelectricPoint(seq_aa, "DTASelect"), 9.645, 
                 tolerance = 1e-03)
    expect_equal(calculateIsoelectricPoint(seq_aa, "Solomon"), 10.046, 
                 tolerance = 1e-03)
    expect_equal(calculateIsoelectricPoint(seq_aa, "Sillero"), 9.921, 
                 tolerance = 1e-03)
    expect_equal(calculateIsoelectricPoint(seq_aa, "Rodwell"), 10.498, 
                 tolerance = 1e-03)
    expect_equal(calculateIsoelectricPoint(seq_aa, "Lehninger"), 10.017, 
                 tolerance = 1e-03)
    expect_equal(calculateIsoelectricPoint(seq_aa, "Toseland"), 9.824,
                 tolerance = 1e-03)
    expect_equal(calculateIsoelectricPoint(seq_aa, "Thurlkill"), 9.866, 
                 tolerance = 1e-03)
    expect_equal(calculateIsoelectricPoint(seq_aa, "Nozaki"), 9.784,
                 tolerance = 1e-03)
    expect_equal(calculateIsoelectricPoint(seq_aa, "IPC_protein"), 8.958, 
                 tolerance = 1e-03)
    expect_equal(calculateIsoelectricPoint(seq_aa, "IPC_peptide"), 10.046, 
                 tolerance = 1e-03)
    
    ## Q9BVP2-1 experimental isoelectric point 8.31
    seq_aa <- AAString("MKRPKLKKASKRMTCHKRYKIQKKVREHHRKLRKEAKKRGHKKPRKDPGVPNSAPFKEALLREAELRKQRLEELKQQQKLDRQKELEKKRKLETNPDIKPSNVEPMEKEFGLCKTENKAKSGKQNSKKLYCQELKKVIEASDVVLEVLDARDPLGCRCPQVEEAIVQSGQKKLVLILNKSDLVPKENLESWLNYLKKELPTVVFRASTKPKDKGKITKRVKAKKNAAPFRSEVCFGKEGLWKLLGGFQETCSKAIRVGVIGFPNVGKSSIINSLKQEQMCNVGVSMGLTRSMQVVPLDKQITIIDSPSFIVSPLNSSSALALRSPASIEVVKPMEAASAILSQADARQVVLKYTVPGYRNSLEFFTVLAQRRGMHQKGGIPNVEGAAKLLWSEWTGASLAYYCHPPTSWTPPPYFNESIVVDMKSGFNLEELEKNNAQSIRAIKGPHLANSILFQSSGLTNGIIEEKDIHEELPKRKERKQEEREDDKDSDQETVDEEVDENSSGMFAAEETGEALSEETTAGEQSTRSFILDKIIEEDDAYDFSTDYV")
    expect_equal(calculateIsoelectricPoint(seq_aa, "EMBOSS"), 9.746,
                 tolerance = 1e-03)
    expect_equal(calculateIsoelectricPoint(seq_aa, "DTASelect"), 9.158, 
                 tolerance = 1e-03)
    expect_equal(calculateIsoelectricPoint(seq_aa, "Solomon"), 9.551, 
                 tolerance = 1e-03)
    expect_equal(calculateIsoelectricPoint(seq_aa, "Sillero"), 9.509,
                 tolerance = 1e-03)
    expect_equal(calculateIsoelectricPoint(seq_aa, "Rodwell"), 10.048, 
                 tolerance = 1e-03)
    expect_equal(calculateIsoelectricPoint(seq_aa, "Lehninger"), 9.533, 
                 tolerance = 1e-03)
    expect_equal(calculateIsoelectricPoint(seq_aa, "Toseland"), 9.33, 
                 tolerance = 1e-03)
    expect_equal(calculateIsoelectricPoint(seq_aa, "Thurlkill"), 9.416,
                 tolerance = 1e-03)
    expect_equal(calculateIsoelectricPoint(seq_aa, "Nozaki"), 9.518, 
                 tolerance = 1e-03)
    expect_equal(calculateIsoelectricPoint(seq_aa, "IPC_protein"), 8.306, 
                 tolerance = 1e-03)
    expect_equal(calculateIsoelectricPoint(seq_aa, "IPC_peptide"), 9.555, 
                 tolerance = 1e-03)
})
