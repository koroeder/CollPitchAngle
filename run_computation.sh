for ((i=1;i<=6275;i++))
do
echo "Structure ${i}" >> output
python ../../TensorComputations/get_tensor.py model_${i}.pdb >> output
rm model_${i}.pdb
done
