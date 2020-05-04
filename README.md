# Code

for specific data downloads from ncbi: https://pypi.org/project/ncbi-genome-download/

# Scripts_short_codes

- exonsInAGene.py gtf file a gore genlerin icindeki exonlari bulmak icin yazildi. Sonrasinda gtfToBed calistirilirsa coordinate leri cekebiliyoruz.

- gtfToBed.py exonlari verilen bir gene in coordinates lerini bulmak icin yazildi. Gene name e gore siralandi.

- gtGeneMapping.py belirli coordinate lerdeki genotype datalarinin hangi gene lere ait oldugunu bulmak icin yazildi.

# Shell commands

* how to login your remote server without password:* >> sudo apt-get install sshpass; alias name_of_command_shortcut='sshpass -p ssh username@remote_host'

* how to change command name:* >> nano ~/.bashrc; alias my_command='original command'; source ~/.bashrc

* how to reach a tool from anywhere in terminal:*  >> nano ~/.bashrc; export PATH=$PATH:/path/; source ~/.bashrc

* install bedtools:* sudo apt install bedtools
