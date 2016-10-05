# This commands should result in an error because the file test01 is empty
python clip/clip.py extract -i tests/testExtract/test01.bam  -c s-1i
python clip/clip.py extract -i tests/testExtract/test01.bam  -c s1i
python clip/clip.py extract -i tests/testExtract/test01.bam  -c s
python clip/clip.py extract -i tests/testExtract/test01.bam  -c e
python clip/clip.py extract -i tests/testExtract/test01.bam  -c e-1i
python clip/clip.py extract -i tests/testExtract/test01.bam  -c e+1i
python clip/clip.py extract -i tests/testExtract/test01.bam  -c m
python clip/clip.py extract -i tests/testExtract/test01.bam  -c i
python clip/clip.py extract -i tests/testExtract/test01.bam  -c d

# This command should result in an error because "bla" is not specified
python clip/clip.py extract -i tests/testExtract/test01.bam  -c bla

# This commands should result in an error because the file test01 is empty
python clip/clip.py extract -i tests/testExtract/test02.bam  -c e

python clip/clip.py extract -i tests/testExtract/test02.bam  -c s-1i
python clip/clip.py extract -i tests/testExtract/test02.bam  -c s1i
python clip/clip.py extract -i tests/testExtract/test02.bam  -c s

python clip/clip.py extract -i tests/testExtract/test02.bam  -c e-1i
python clip/clip.py extract -i tests/testExtract/test02.bam  -c e+1i
python clip/clip.py extract -i tests/testExtract/test02.bam  -c m
python clip/clip.py extract -i tests/testExtract/test02.bam  -c i
python clip/clip.py extract -i tests/testExtract/test02.bam  -c d