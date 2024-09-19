
source ~/.bashrc
conda activate intervene



intervene upset -i vRG.Both.peak.bed oRG.Both.peak.bed OPC.Both.peak.bed MG.Both.peak.bed \
                    --output ./Intervene_Both --save-overlaps \
                    --names vRG,oRG,OPC,MG \
                    --figsize 5 3
                    
                    
intervene venn -i vRG.Both.peak.bed oRG.Both.peak.bed OPC.Both.peak.bed MG.Both.peak.bed \
                    --output ./Intervene_Both --save-overlaps \
                    --names vRG,oRG,OPC,MG \
                    --colors='#2166AC','#B2182B','#1B7837','#762A83' \
                    --figsize 12 5 \
                    --fontsize 7