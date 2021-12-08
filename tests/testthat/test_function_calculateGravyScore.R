## sample sequences taken from
## https://www.bioinformatics.org/sms2/protein_gravy.html
test_that("calculateGravyScore", {
    
    expect_error(calculateGravyScore(1),
                 "'aa' has to be a character vector")
    expect_error(calculateGravyScore(c("MQ", "AELS")),
                 "'aa' has to be of length 1")
    expect_equal(calculateGravyScore("xycAME"), NaN)
    seq <- "MQKSPLEKASFISKLFFSWTTPILRKGYRHHLELSDIYQAPSADSADHLSEKLEREWDREQASKKNPQLIHALRRCFFWRFLFYGILLYLGEVTKAVQPVLLGRIIASYDPENKVERSIAIYLGIGLCLLFIVRTLLLHPAIFGLHRIGMQMRTAMFSLIYKKTLKLSSRVLDKISIGQLVSLLSNNLNKFDEGLALAHFIWIAPLQVTLLMGLLWDLLQFSAFCGLGLLIILVIFQAILGKMMVKYRDQRAAKINERLVIT"
    expect_equal(calculateGravyScore(seq), 0.308, tolerance = 1e-03)
    seq <- "SEIIDNIYSVKAYCWESAMEKMIENLREVELKMTRKAAYMRFFTSSAFFFSGFFVVFLSVLPYTVINGIVLRKIFTTISFCIVLRMSVTRQFPTAVQIWYDSFGMIRKIQDFLQKQEYKVLEYNLMTTGI"
    expect_equal(calculateGravyScore(seq), 0.321, tolerance = 1e-03)
    seq <- "IMENVTAFWEEGFGELLQKAQQSNGDRKHSSDENNVSFSHLCLVGNPVLKNINLNIEKGEMLAITGSTGLGKTSLLMLILGELEASEGIIKHSGRVSFCSQFSWIMPGTIKENIIFGVSYDEYRYKSV"
    expect_equal(calculateGravyScore(seq), -0.127, tolerance = 1e-03)
})