# Improving gRNA efficiency of CRISPR/Cpf1 system
This study was conducted to develop a off-target-prediction pipeline and study biological importance of features using machine learning approach and statistical analyses. Multiple classification algorithms were trained on two different datasets (AsCpf1 and LbCpf1) and their performances was compared. Best performing models was used in the off-target prediction pipelines. 

## Feature importance analysis
Statisitical significance of each position-specific nucleotides were estimated to understand their biological importance in the classification of off-targets. Other mismatch associated features were also explored and their importances were studied using statistical analysis. Machine learning approach was also used to analyze the importance of features. All the neccessary analysis codes are given in Jupyter Notebooks.

## off-target prediction pipeline
Off-target prediction pipeline is divided into two parts which are to calculate the target efficiency and predict potential off-target sites in Human genome. 

### Usage
To implement the off-target prediction pipeline in your linux-based system: <br />
Reqirements: <br />
  •	FASTA36 <br />
  •	ViennaRNA package <br />
  •	SeqIO, SeqUtils, Seq modules of Bio package <br />
  •	Biopython <br />
  •	Python packages: Pandas, Numpy, sci-kit learn, joblib <br />
  <br />
Execution <br />
 • For Potential off-target site prediction: <br />
  <br />
  Python3 Fasta36_SeqSearch.py input_filename.fasta <br />
  <br />
  • For target efficiency prediction: <br />
  <br />
  Python3 main.py input_filename.fasta <br />
<br />
Refer usage file for detail implementation guide for the off-target prediction pipeline

### Have suggestions or need help?
In case of any queries and suggestions, please [submit a query or an issue.](https://github.com/PragyaKesar/CRISPR-Cpf1_study/issues/new)

(Suggestions are always welcome! :smiley: )
