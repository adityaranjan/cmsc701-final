python run_all.py LexMin 1 0
mv data/results/* results/lexmin/partial_1
python run_all.py LexMin 2 0
mv data/results/* results/lexmin/partial_2
python run_all.py LexMin 3 0
mv data/results/* results/lexmin/partial_3
python run_all.py LexMin 4 0
mv data/results/* results/lexmin/partial_4

python run_all.py LexMin 0 1
mv data/results/* results/lexmin/delta_1
python run_all.py LexMin 0 2
mv data/results/* results/lexmin/delta_2
python run_all.py LexMin 0 3
mv data/results/* results/lexmin/delta_3
python run_all.py LexMin 0 4
mv data/results/* results/lexmin/delta_4

python run_all.py LexMin 1 1
mv data/results/* results/lexmin/combo_1_1
python run_all.py LexMin 1 2
mv data/results/* results/lexmin/combo_1_2
python run_all.py LexMin 2 1
mv data/results/* results/lexmin/combo_2_1
python run_all.py LexMin 2 2
mv data/results/* results/lexmin/combo_2_2

python run_all.py HashMin 1 0
mv data/results/* results/hashmin/partial_1
python run_all.py HashMin 2 0
mv data/results/* results/hashmin/partial_2
python run_all.py HashMin 3 0
mv data/results/* results/hashmin/partial_3
python run_all.py HashMin 4 0
mv data/results/* results/hashmin/partial_4

python run_all.py HashMin 0 1
mv data/results/* results/hashmin/delta_1
python run_all.py HashMin 0 2
mv data/results/* results/hashmin/delta_2
python run_all.py HashMin 0 3
mv data/results/* results/hashmin/delta_3
python run_all.py HashMin 0 4
mv data/results/* results/hashmin/delta_4

python run_all.py HashMin 1 1
mv data/results/* results/hashmin/combo_1_1
python run_all.py HashMin 1 2
mv data/results/* results/hashmin/combo_1_2
python run_all.py HashMin 2 1
mv data/results/* results/hashmin/combo_2_1
python run_all.py HashMin 2 2
mv data/results/* results/hashmin/combo_2_2