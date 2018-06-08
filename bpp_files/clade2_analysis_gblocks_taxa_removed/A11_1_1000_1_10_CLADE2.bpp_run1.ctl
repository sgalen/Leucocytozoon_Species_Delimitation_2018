          seed =  1

       seqfile = /home/sgalen/clade2_analysis_gblocks/clade2_gblocks.phy
      Imapfile = /home/sgalen/clade2_analysis_gblocks/leuco_imap_clade2.txt
       outfile = /home/sgalen/clade2_analysis_gblocks/A11_1_1000_1_10_clade2_out_run1.txt
      mcmcfile = /home/sgalen/clade2_analysis_gblocks/A11_1_1000_1_10_clade2_mcmc_run1.txt

*  speciesdelimitation = 0 * fixed species tree
 speciesdelimitation = 1 0 2   * species delimitation rjMCMC algorithm0 and finetune(e)
* speciesdelimitation = 1 1 2 1  * species delimitation rjMCMC algorithm1 finetune (a m)
         speciestree = 1  *0.4 0.2 0.1   * speciestree pSlider ExpandRatio ShrinkRatio

   speciesmodelprior = 1  * 0: uniform LH; 1:uniform rooted trees; 2: uniformSLH; 3: uniformSRooted

  species&tree = 11  CATMIN01  PERCAN01  ACAFLA03  SETPET04  CATUST14  BOMGAR01  ROFI06  LOXLEU05  CNEORN01  ZOLEU02  ZOLEU04
                    1  6  6  1  2  1  3  1  5  6  1  
                   (CATMIN01,(((PERCAN01,((ACAFLA03,SETPET04),CATUST14)),(BOMGAR01,(ROFI06,LOXLEU05))),(CNEORN01,(ZOLEU02,ZOLEU04))));




                  
       usedata = 1  * 0: no data (prior); 1:seq like
         nloci = 7  * number of data sets in seqfile

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

    thetaprior = 1 1000   # gamma(a, b) for theta
      tauprior = 1 10   # gamma(a, b) for root tau & Dirichlet(a) for other tau's

*     heredity = 1 4 4
*    locusrate = 1 5

      finetune =  1: 5 0.001 0.001  0.001 0.3 0.33 1.0  # finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr

         print = 1 0 0 0   * MCMC samples, locusrate, heredityscalars, Genetrees
        burnin = 25000
      sampfreq = 5
       nsample = 250000
