

# Test:

## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## Section: headerReadGroup
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## Create headerReadGroup by simply parsing text
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

hrg <- c("@RG\tID:rg1\tCN:seqCenter1\tDT:01.01.2011\tSM:sm1",
                    "@RG\tID:rg2\tCN:seqCenter2\tDT:01.01.2012\tSM:sm2")
object <- new("headerReadGroup", hrg)
identical(object@ID,paste("rg", 1:2, sep=""))


## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## Create headerReadGroup by user interface
## and check for consistency of returned header text
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

hrg <- new("headerReadGroup")
hrg <- addReadGroup(hrg, list(ID="rg1", CN="sct1", FO="fo1"))
hrg <- addReadGroup(hrg, list(ID="rg2", CN="sct2", FO="fo2", LB="lb2"))
hrg <- addReadGroup(hrg, list(ID="rg3", CN="sct3", LB="lb3"))
hrg2 <- new("headerReadGroup", getHeaderText(hrg))

if(!identical(hrg, hrg2))
    stop("[test_bam_header.r] headerReadGroup test failed")
