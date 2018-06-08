          seed =  1

       seqfile = /home/sgalen/clade1_catust_tumig_removed_gblocks/clade1_catust09.tumig15_removed.phy
      Imapfile = /home/sgalen/clade1_catust_tumig_removed_gblocks/leuco_imap_clade1_catust09_tumig15_removed.txt
       outfile = /home/sgalen/clade1_catust_tumig_removed_gblocks/A11_1_1000_1_1000_clade1_catust09.tumig15_removed_out_run1.txt
      mcmcfile = /home/sgalen/clade1_catust_tumig_removed_gblocks/A11_1_1000_1_1000_clade1_catust09.tumig15_removed_mcmc_run1.txt

*  speciesdelimitation = 0 * fixed species tree
 speciesdelimitation = 1 0 2   * species delimitation rjMCMC algorithm0 and finetune(e)
* speciesdelimitation = 1 1 2 1  * species delimitation rjMCMC algorithm1 finetune (a m)
         speciestree = 1  *0.4 0.2 0.1   * speciestree pSlider ExpandRatio ShrinkRatio

   speciesmodelprior = 1  * 0: uniform LH; 1:uniform rooted trees; 2: uniformSLH; 3: uniformSRooted

  species&tree = 15  JUNHYE04  CARPUS01  CB1  ACAFLA04  PHYBOR01  BT1  ACAFLA01  TRPIP2  PERCAN04  COBRA07  TUMIG20  TUMIG11   CATMIN05  CATUST28  CATGUT02  
                    1  1  10  1  1  3  1  1  1  1  1  1  1  1  5  
(((((((JUNHYE04,CARPUS01),CB1),ACAFLA04),(PHYBOR01,BT1)),(ACAFLA01,TRPIP2)),(PERCAN04,COBRA07)),((TUMIG20,TUMIG11),((CATMIN05,CATUST28),CATGUT02)));




                  
       usedata = 1  * 0: no data (prior); 1:seq like
         nloci = 6  * number of data sets in seqfile

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

    thetaprior = 1 1000   # gamma(a, b) for theta
      tauprior = 1 1000   # gamma(a, b) for root tau & Dirichlet(a) for other tau's

*     heredity = 1 4 4
*    locusrate = 1 5

      finetune =  1: 5 0.001 0.001  0.001 0.3 0.33 1.0  # finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr

         print = 1 0 0 0   * MCMC samples, locusrate, heredityscalars, Genetrees
        burnin = 25000
      sampfreq = 5
       nsample = 250000
