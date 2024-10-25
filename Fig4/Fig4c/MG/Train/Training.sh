#bin/bash
#$ -l h_rt=24:0:0
#$ -l mem_free=16G
#$  -cwd
#$ -l scratch=500G
#$ -j y
#$ -r y

#train the regression model
#-y 3 selects regression (default of -y 0 does classification)
#-t controls the choice of kernel. -t 2 (instead of the default of 4)
# avoids upweighting the central positions.
#The parameter -p controls the “epsilon” used in the support vector
# regression objective (i.e. the tolerated margin of error).
# HEADS UP: don’t confuse this with the -e parameter, which determines
# the epsilon used for the convergence criterion (they are both
# called ‘epsilon’, sigh).
#Depending on the problem, the -c parameter (controlling the cost of
# misprediction) will likely have to be adjusted, i.e. the defaults may
# not work out-of-the-box.
/shen/shenlabstore3/ijones1/dependencies/GKM_Explain/rlsgkm/bin/gkmtrain -r 1 -T 16 -y 0 -t 4 -c 1 -p 0.1 /shen/shenlabstore3/ijones1/GKM_explain_test/Microglia/Pos_90k/MG.train.75k.peaks.filt.fa /shen/shenlabstore3/ijones1/GKM_explain_test/Microglia/Neg_90k/MG.train.75k.neg.filt.fa MG.Not_chr2.model.8.5.24