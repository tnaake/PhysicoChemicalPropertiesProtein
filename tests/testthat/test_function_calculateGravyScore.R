## sample sequences taken from
## https://www.bioinformatics.org/sms2/protein_gravy.html
test_that("calculateGravyScore works.", {
    
    expect_error(calculateGravyScore(1),
                 "'aa' has to be an 'AAString' object")
    expect_error(calculateGravyScore(c("MQ", "AELS")),
                 "'aa' has to be an 'AAString' object")
    seq_aa <- AAString("xycAME")
    expect_error(calculateGravyScore(seq_aa), 
                 "'aa' contains letters not defined in the amino acid alphabet")
    seq_aa <- AAString("MQKSPLEKASFISKLFFSWTTPILRKGYRHHLELSDIYQAPSADSADHLSEKLEREWDREQASKKNPQLIHALRRCFFWRFLFYGILLYLGEVTKAVQPVLLGRIIASYDPENKVERSIAIYLGIGLCLLFIVRTLLLHPAIFGLHRIGMQMRTAMFSLIYKKTLKLSSRVLDKISIGQLVSLLSNNLNKFDEGLALAHFIWIAPLQVTLLMGLLWDLLQFSAFCGLGLLIILVIFQAILGKMMVKYRDQRAAKINERLVIT")
    expect_equal(calculateGravyScore(seq_aa), 0.308, tolerance = 1e-03)
    seq_aa <- AAString("SEIIDNIYSVKAYCWESAMEKMIENLREVELKMTRKAAYMRFFTSSAFFFSGFFVVFLSVLPYTVINGIVLRKIFTTISFCIVLRMSVTRQFPTAVQIWYDSFGMIRKIQDFLQKQEYKVLEYNLMTTGI")
    expect_equal(calculateGravyScore(seq_aa), 0.321, tolerance = 1e-03)
    seq_aa <- AAString("IMENVTAFWEEGFGELLQKAQQSNGDRKHSSDENNVSFSHLCLVGNPVLKNINLNIEKGEMLAITGSTGLGKTSLLMLILGELEASEGIIKHSGRVSFCSQFSWIMPGTIKENIIFGVSYDEYRYKSV")
    expect_equal(calculateGravyScore(seq_aa), -0.127, tolerance = 1e-03)
})
