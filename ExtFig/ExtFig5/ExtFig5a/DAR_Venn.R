dsource ~/.bashrc
conda activate intervene


intervene venn -i /shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/Ziffra_data/tRG.GT.oRG.Ziffra.bed /shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/Ziffra_data/oRG.GT.tRG.Ziffra.bed vRG.DAR.bed oRG.DAR.bed \
                    --output ./Intervene_Zifra --save-overlaps \
                    --names tRG_Ziffra,oRG_Ziffra,vRG,oRG \
                    --colors='#92C5DE','#F4A582','#2166AC','#B2182B' \
                    --figsize 6 4 \
                    --fontsize 7