# Trabalho final - Pipeline Somático 2
Trabalho final do curso de pós-graduação em Bioinformática Aplicada a Genômica Médica - Einstein

Integrantes:
- Ana Vitória Volpato Jensen
- Angelica Nakagawa Lima
- Antônio José Sousa Soares
- Edgar da Costa Silva
- Nicholy Cristine Forner Lozano

## Filtragem de variantes somáticas - Painel de Mielofibrose com prognóstico ruim
 - Este pipeline é destinado a filtrar variantes somáticas que ofereçam prognóstico ruim para a doença de mielofibrose.
 - Foram filtradas as amostras do projeto "LMA Brasil", utilizando os seguintes genes para a realizar a filtragem:
*CALR*, *CBL*, *EZH2*, *GFI1B*, *IDH1*, *IDH2*, *JAK2*, *KRAS*, *MPIG6B*, *MPL*, *NBEAL2*, *NRAS*, *SH2B3*, *SHOC2*, *SRC*, *SRSF2*, *TBXAS1*, *TET2*, *TLR8*, *TP53*, *U2AF1*.

![image](https://github.com/angelicanakagawa/TrabalhoFinal_Somatico2/assets/91493865/6d05c667-7227-4f89-9d2a-fe8061b36701)
**Figura 1.** Esquema do pipeline realizado neste trabalho. Autoria: Ana Vitória V. Jensen

### Etapa 1. Preparar do ambiente

1.1. Obter acesso aos arquivos dentro de seu Drive através do Colab

```
from google.colab import drive
drive.mount('/content/drive')
```

1.2. Clonar o github do projeto lmabrasil

```
%%bash
rm -rf lmabrasil-hg38
git clone https://github.com/renatopuga/lmabrasil-hg38
```

1.3. Clonar o github do projeto ***TrabalhoFinal_Somatico2***

```
%%bash
git clone https://github.com/angelicanakagawa/TrabalhoFinal_Somatico2
```

1.4. Instalar o BCFTools + split-vep

```
%%bash
git clone --recurse-submodules https://github.com/samtools/htslib.git
git clone https://github.com/samtools/bcftools.git
cd bcftools
make
make install
```

1.5. Instalar o GATK - download e descompactação

```
%%bash
wget -c https://github.com/broadinstitute/gatk/releases/download/4.2.6.1/gatk-4.2.6.1.zip
unzip gatk-4.2.6.1.zip
```

1.6. Instalar o tabix

```
%%bash
apt-get install tabix
```

1.7. Instalar o UDocker. 
> Fonte: https://gist.github.com/mwufi/6718b30761cd109f9aff04c5144eb885; https://github.com/indigo-dc/udocker

```
%%bash
pip install udocker
udocker --allow-root install
```

1.8. Download da imagem do Ensembl-VEP utilizada no UDocker

```
%%bash
udocker --allow-root pull ensemblorg/ensembl-vep
```

### Etapa 2. Realizar o LiftOver

Os arquivos deste projeto estão em hg19 (GRCh37), sendo necessário transformá-los em hg38 (GRCh38), através do liftover

2.1. Download dos arquivos VCFs na versão hg19

> Fonte: https://drive.google.com/drive/folders/1m2qmd0ca2Nwb7qcK58ER0zC8-1_9uAiE?usp=sharing

```
%%bash
wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
```

2.2. Descompactar o arquivo baixado

```
%%bash
gunzip hg19ToHg38.over.chain.gz
```

2.3. Colocar as nomenclaturas dos cromossomos no mesmo formato ("chr1, chr2, chr3, ...")

> Os VCF's estão escritos com diferentes nomenclaturas para a sinalização dos cromossomos. Nos arquivos do projeto ***TrabalhoFinal_Somatico2*** estão no formato "1, 2, 3, ...", enquanto que no arquivo que será utilizado no LiftOver está como "chr1, chr2, chr3, ..."

```
%%bash
# Criar o caminho dos VCF's deste projeto
path_vcf="/content/TrabalhoFinal_Somatico2/arquivos_lmabrasil/vcfs_hg19"

# Criar diretório de resultados
mkdir resultados
mkdir resultados/vcfs_hg38

# Entrar no diretório de resultados
cd resultados/vcfs_hg38

# Loop para alterar todos os vcfs
for filename in ${path_vcf}/*.filtered.vcf.gz; do
  # Extrair o nome da amostra
  sample=$(basename "$filename" .vcf.gz)

  # Extrair o cabeçalho
  zgrep "\#" ${path_vcf}/${sample}.vcf.gz > header.txt

  # Adicionar o "chr" na primeira coluna
  zgrep -v "\#" ${path_vcf}/${sample}.vcf.gz | awk '{print("chr"$0)}' > variants.txt

  # Incluir o cabeçalho
  cat header.txt variants.txt > ${sample}.chr.vcf
done
```

2.4. Compactar os arquivos VCF's com bgzip e criando os índices

```
%%bash
# Entrar no diretório de resultados
cd resultados/vcfs_hg38

# Loop para executar todos os vcfs:
for filename in *.filtered.chr.vcf; do

   # Compactar o vcf
   bgzip ${filename}

   # Criar o índice
   tabix -p vcf ${filename}.gz
done
```

2.5. Baixar o genoma de referência hg38

```
%%bash
# Baixar os arquivos
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dict

# Criar um diretório para armazenar os genomas de referência
mkdir genoma_referencia

# Mover os arqivos para o diretório de genomas de referência
mv Homo_sapiens_assembly38* genoma_referencia/
```

2.6. Fazendo o LiftoverVcf usando o gatk

```
%%bash
# Entrar no diretório de resultados
cd resultados/vcfs_hg38

# Loop para executar todos os VCF's:
for filename in *.filtered.chr.vcf.gz; do

  # Extrair os nomes das amostras
  sample = $(basename "$filename" .filtered.chr.vcf.gz)

  # Realizar o LiftOver dos VCF's através do GATK
  /content/gatk-4.2.6.1/gatk LiftoverVcf \
   -I ${sample}.filtered.chr.vcf.gz \
   -O liftOver_${sample}\_hg19ToHg38.vcf \
   --CHAIN /content/hg19ToHg38.over.chain \
   --REJECT liftOver_Reject_${sample}\_hg19ToHg38.vcf \
   -R /content/genoma_referencia/Homo_sapiens_assembly38.fasta
done
```

### Etapa 3. Anotar os VCF's

Esta etapa é demorada para ser realizada no colab. Portanto, os comandos foram executados localmente e os VCF's anotados foram copiados para o Drive.

Abaixo estão os comandos utilizados para uma única amostra (WP306):

```
%%bash
# Criar um diretório para armazenar os resultados finais da anotação (utilizar o comando "chmod 777" para permitir que usuários leiam, editem e executem arquivos no diretório criado)
mkdir -p vep_output
chmod 777 vep_output

# Baixar a referência para ser usada no VEP
curl -O https://ftp.ensembl.org/pub/release-110/variation/indexed_vep_cache/homo_sapiens_vep_110_GRCh38.tar.gz
tar xzf homo_sapiens_vep_110_GRCh38.tar.gz
cd resultados
REF="/content/genoma_referencia"

# Realizar o LiftOver
udocker run --allow-root  -it --rm -v $(pwd):/data -v $REF:/referencia ensemblorg/ensembl-vep:release_110.1 vep \
-i /data/liftOver_WP306_hg19ToHg38.vcf  \
-o /data/vep_output/liftOver_WP306_hg19ToHg38.vcf \
--assembly GRCh38  \
--refseq \
--fork 8 \
--buffer_size 200 \
--force_overwrite \
--dir_cache /data \
--offline \
--cache \
--no_intergenic \
--distance 0 \
--pick \
--individual all \
--vcf \
--symbol \
--biotype \
--hgvs \
--mane_select \
--numbers \
--af_gnomadg \
--max_af \
--variant_class \
--sift b \
--variant_class \
--polyphen b \
--check_existing \
--fields "Location,SYMBOL,Consequence,Feature,MANE_SELECT,BIOTYPE,HGVSc,HGVSp,EXON,INTRON,VARIANT_CLASS,SIFT,PolyPhen,gnomADg_AF,MAX_AF,IMPACT,CLIN_SIG,SOMATIC,Existing_variation" \
--fasta /REF/Homo_sapiens_assembly38.fasta
```

### Etapa 4. Filtragem

4.1. Fazer uma cópia dos VCF's já anotados, contendo as amostras pós LiftOver para hg38

```
%%bash
# Remover o conteúdo existente na pasta "vep_output" (caso tenha sido feito o teste com a amostra "WP306")
rm /content/lmabrasil-hg38/vep_output/*

# Copiar os 30 arquivos da pasta "projetolma" (contendo amostras pós lift-over hg38) para a pasta "vep_output"
cp /content/TrabalhoFinal_Somatico2/arquivos_lmabrasil/lmabrasil-lifOverhg19ToHg38/* /content/lmabrasil-hg38/vep_output
```

4.2. Criar uma lista com 9 genes de risco (+ 12 genes adicionais) para prognóstico ruim para mielofibrose

```
%%bash
echo -e "TP53\nCALR\nGFI1B\nJAK2\nMPIG6B\nMPL\nNBEAL2\nSH2B3\nSHOC2\nSRC\nTBXAS1\nTET2\nTLR8\nEZH2\nCBL\nU2AF1\nSRSF2\nIDH1\nIDH2\nNRAS\nKRAS\n" > /content/lmabrasil-hg38/hpo/mielofibrose.txt
```

4.3. Listar os nomes das 30 amostras do projeto LMA Brasil

```
SAMPLES = ["WP048","WP093","WP087","WP060","WP056","WP066","WP064","WP072","WP078","WP285","WP280","WP274","WP276","WP270","WP216","WP306","WP297","WP291","WP295","WP204","WP160","WP164","WP162","WP212","WP170","WP196","WP180","WP188","WP140","WP126"]
```

4.4. Descompactação dos arquivos VCF's

```
%%bash

# Entrar no diretório de resultados
cd lmabrasil-hg38/vep_output/

# Loop para executar todos os vcfs:
for filename in *.vep.vcf.gz; do
  # extraindo o nome da amostra
  sample=$(basename "$filename" .vep.vcf.gz)
  # descompactação
  gunzip $sample.vep.vcf.gz
  # remove arquivos .tbi
  rm $sample.vep.vcf.gz.tbi
done
```

4.5. Loop para processar pipeline de filtragem de variantes em cada uma das 30 amostras da lista SAMPLES

```
for i in SAMPLES:
  SAMPLE = i
  !echo {SAMPLE}
  !sh lmabrasil-hg38/vep-gc.sh {SAMPLE} mielofibrose.txt
```

4.6. Converter as 30 amostras filtradas pelo VEP ("vep-gc.sh") em amostras ".csv" e coloca-lás em um output próprio

```
import pandas as pd
pd.set_option('display.max_columns', None)
pd.set_option('display.max_colwidth', None)

# Criar um diretório para as amostras
!mkdir /content/lmabrasil-hg38/csv_filtrados

# Colocar as amostras no local certo e em ".csv"
for i in SAMPLES:
  df = pd.read_csv(f'/content/lmabrasil-hg38/vep_output/liftOver_{i}_hg19ToHg38.vep.filter.tsv',sep='\t',index_col=False)
  df.to_csv(f'/content/lmabrasil-hg38/csv_filtrados/{i}_filtrado.csv', index = False)
```

### Etapa 5. Montar uma tabela final e vizualizar os resultados

Resultado da filtragem das 30 amostras para mutações de prognóstico ruim para mielofibrose:
  - **Total de variantes encontradas nod genes alvo: 9**
  - **Amostras com variantes encontradas nos genes alvo (n = 8): WP164, WP060, WP216, WP048, WP212, WP306, WP270, WP276**

```
import glob
import pandas as pd

# Listar pasta com os 30 arquivos CSV's
csv_files = glob.glob('/content/lmabrasil-hg38/csv_filtrados/*.{}'.format('csv'))

# Unir o resultado dos 30 arquivos CSV's numa única tabela
df_concat = pd.concat([pd.read_csv(i) for i in csv_files], ignore_index=True)

# Criar um diretório para a tabela final transformando a mesma para o formato CSV
!mkdir /content/lmabrasil-hg38/tabela_final
df_concat.to_csv('/content/lmabrasil-hg38/tabela_final/tabela_final.csv', index = False)
df_concat

# Alterar a coluna "TumorId" para "AMOSTRA"
dados = pd.read_csv('/content/lmabrasil-hg38/tabela_final/tabela_final.csv')
dados = dados.rename(columns={'TumorID': 'AMOSTRA'})
dados.to_csv('/content/lmabrasil-hg38/tabela_final/tabela_final.csv', index = False)
```

Figura 1 - Quantidade de variantes por gene

```
import pandas as pd
import matplotlib.pyplot as plt
import plotly.express as px

# Criar o caminho da tabela final
dados = pd.read_csv('/content/lmabrasil-hg38/tabela_final/tabela_final.csv')

# Agrupar os dados por 'AMOSTRA' e 'SYMBOL' e contar o número de ocorrências de cada combinação
dados_agrupados = dados.groupby(['AMOSTRA', 'SYMBOL']).size().unstack(fill_value=0)

# Resetar o índice para tornar 'AMOSTRA' uma coluna
dados_agrupados = dados_agrupados.reset_index()

# Reorganizar o DataFrame para usa-lo no Plotly Express
dados_melted = pd.melt(dados_agrupados, id_vars = 'AMOSTRA', var_name = 'Variante', value_name = 'Quantidade')

# Criar um gráfico de barras através do Plotly Express
fig = px.bar(dados_melted, x = 'AMOSTRA', y = 'Quantidade', color = 'Variante', barmode = 'stack', title = 'Variantes por Amostra')
fig.update_layout(xaxis_title = 'Amostras', yaxis_title = 'Quantidade de Variantes')
#exporta o Grafico pra html
fig.write_html('/content/resultados/Variantes.html')

fig.show()
```

Figura 2 - Variantes de alto risco

```
import pandas as pd
import matplotlib.pyplot as plt
import plotly.express as px

# Criar o caminho da tabela final
dados = pd.read_csv('/content/lmabrasil-hg38/tabela_final/tabela_final.csv')

# Filtrar os genes de alto risco
genes_filtrados = dados[dados['SYMBOL'].isin(['KRAS', 'SRSF2', 'U2AF1', 'NRAS', 'TP53'])]

# Agrupar os dados por 'AMOSTRA' e 'SYMBOL' e contar o número de ocorrências de cada combinação
dados_agrupados = genes_filtrados.groupby(['AMOSTRA', 'SYMBOL']).size().unstack(fill_value=0)

# Resetar o índice para tornar 'AMOSTRA' uma coluna
dados_agrupados = dados_agrupados.reset_index()

# Melt do DataFrame para o formato necessário para Plotly Express
dados_melted = pd.melt(dados_agrupados, id_vars = 'AMOSTRA', var_name = 'Variante', value_name = 'Quantidade')

# Criar gráfico de barras interativo com Plotly Express
fig = px.bar(dados_melted, x = 'AMOSTRA', y = 'Quantidade', color = 'Variante', barmode = 'stack', title = 'Genes de alto risco', hover_name = 'AMOSTRA')
fig.update_layout(xaxis_title = 'Amostra/Paciente', yaxis_title = 'Quantidade de Variantes')
fig.write_html('/content/resultados/VariantesAltoRisco.html')
```

Figura 3 - Variantes de alto risco

```
import pandas as pd
import matplotlib.pyplot as plt
import plotly.express as px

# Criar o caminho da tabela final
dados = pd.read_csv('/content/lmabrasil-hg38/tabela_final/tabela_final.csv')

# Filtrar os genes de alto risco
genes_filtrados = dados[dados['SYMBOL'].isin(['KRAS', 'SRSF2', 'U2AF1', 'NRAS', 'TP53'])]

# Contar a frequência de cada variantes
contagem_variantes = genes_filtrados['SYMBOL'].value_counts().reset_index()

contagem_variantes.columns = ['Variante', 'Quantidade']

# Criar um gráfico de pizza interativo com Plotly
fig = px.pie(contagem_variantes, values = 'Quantidade', names = 'Variante', title = 'Genes de alto risco')

fig.write_html('/content/resultados/VariantesAltoRiscoPizza.html')
fig.show()
```

Figura 4 - Frequência de amostras com alterações em genes de alto risco

```
import numpy as np
import matplotlib.pyplot as plt

# Frequência de amostras com alterações em genes de alto risco
labels = 'Sem alteração', 'Com alteração'
sections = [21, 9]

plt.pie(sections, labels=labels, autopct = '%1.1f%%')

plt.title('Frequência de amostras com alterações em genes de alto risco')
plt.show()
```

### Encontrar variantes em genes driver

> Curiosamente, dados da literatura (https://doi.org/10.1200/JCO.2016.70.7968) mostram que mais de 90% dos casos de mielofibrose envolvem mutações somáticas nos genes driver *JAK2*, *CALR* ou *MPL*, o que leva a uma ativação constitutiva da via *JAK-STAT5*.

```
import pandas as pd
import matplotlib.pyplot as plt
import plotly.express as px

# Criar o caminho da tabela final
dados = pd.read_csv('/content/lmabrasil-hg38/tabela_final/tabela_final.csv')

# Filtrar os genes drive
genes_filtrados = dados[dados['SYMBOL'].isin(['JAK2', 'CALR','MPL'])]

# Agrupar os dados por 'AMOSTRA' e 'SYMBOL' e contar o número de ocorrências de cada combinação
dados_agrupados = genes_filtrados.groupby(['AMOSTRA', 'SYMBOL']).size().unstack(fill_value = 0)

# Resetar o índice para tornar 'AMOSTRA' uma coluna
dados_agrupados = dados_agrupados.reset_index()

# Melt do DataFrame para o formato necessário para Plotly Express
dados_melted = pd.melt(dados_agrupados, id_vars = 'AMOSTRA', var_name = 'Variante', value_name = 'Quantidade')

# Criar gráfico de barras interativo com Plotly Express
fig = px.bar(dados_melted, x = 'AMOSTRA', y = 'Quantidade', color = 'Variante', barmode = 'stack', title = 'Genes Drive', hover_name = 'AMOSTRA')
fig.update_layout(xaxis_title = 'Amostra/Paciente', yaxis_title = 'Quantidade de Variantes')
fig.write_html('/content/resultados/GenesDrive.html')
fig.show()
```

> Etapa 6. Interpretador do Genoma do Câncer (CGI)

Gerando o arquivos .txt com as colunas: CHR, POS, REF e ALT. Esses arquivos serão utilizandos juntos com a API.

```
%%bash
mkdir /content/CGI

cd /content/lmabrasil-hg38/vep_output/

for file in *.vep.filter.tsv; do
    base=$(basename "$file" .vep.filter.tsv)
    cut -f1-4 "$file" | sed -e "s/CHROM/CHR/g" > "/content/CGI/df_${base}-cgi.txt"
done
```

Enviando um Job
```
import requests
headers = {'Authorization': 'antoniosousa.js98@gmail.com bf6acfa27c682e8b136d'}
payload = {'cancer_type': 'HEMATO', 'title': 'Somatic MF', 'reference': 'hg38'}
r = requests.post('https://www.cancergenomeinterpreter.org/api/v1',
                headers=headers,
                files={
                        'mutations': open('/content/CGI/df_liftOver_WP276_hg19ToHg38-cgi.txt', 'rb')
                        },
                data=payload)
r.json()
```
Visualizando os identificadores
```
import requests
job_id ="149e726abf0bd8ad9f6f"

headers = {'Authorization': 'antoniosousa.js98@gmail.com bf6acfa27c682e8b136d'}
r = requests.get('https://www.cancergenomeinterpreter.org/api/v1/%s' % job_id, headers=headers)
r.json()
```
Acessando informações do Job
```
import requests
job_id ="149e726abf0bd8ad9f6f"

headers = {'Authorization': 'antoniosousa.js98@gmail.com bf6acfa27c682e8b136d'}
payload={'action':'logs'}
r = requests.get('https://www.cancergenomeinterpreter.org/api/v1/%s' % job_id, headers=headers, params=payload)
r.json()
```
Download dos Resultados (file.zip)
```
import requests
job_id ="149e726abf0bd8ad9f6f"

headers = {'Authorization': 'antoniosousa.js98@gmail.com bf6acfa27c682e8b136d'}
payload={'action':'download'}
r = requests.get('https://www.cancergenomeinterpreter.org/api/v1/%s' % job_id, headers=headers, params=payload)
with open('/content/CGI/file.zip', 'wb') as fd:
    fd.write(r._content)
```
```
!unzip /content/CGI/file.zip
```
```
import pandas as pd

pd.read_csv('/content/CGI/alterations.tsv', sep='\t', index_col=False, engine='python')


pd.read_csv('/content/CGI/biomarkers.tsv',sep='\t',index_col=False, engine= 'python')

import requests
#job_id ="3cf2faf653502b3b458d"

headers = {'Authorization': 'antoniosousa.js98@gmail.com bf6acfa27c682e8b136d'}
r = requests.delete('https://www.cancergenomeinterpreter.org/api/v1/%s' % job_id, headers=headers)
r.json()
```

> Teste CGI
```
mkdir /content/CGI
cd /content/lmabrasil-hg38/vep_output/

# Gerar o arquivo "df_WP048-cgi.txt" com as colunas: CHR, POS, REF e ALT
cut -f1-4 /content/lmabrasil-hg38/vep_output/liftOver_WP048_hg19ToHg38.vep.filter.tsv | sed -e "s/CHROM/CHR/g"  > df_WP048-cgi.txt
head df_WP048-cgi.txt

# Enviando um Job
import requests
headers = {'Authorization': 'antoniosousa.js98@gmail.com bf6acfa27c682e8b136d'}
payload = {'cancer_type': 'HEMATO', 'title': 'Somatic MF', 'reference': 'hg38'}
r = requests.post('https://www.cancergenomeinterpreter.org/api/v1',
                headers=headers,
                files={
                        'mutations': open('/content/CGI/df_liftOver_WP276_hg19ToHg38-cgi.txt', 'rb')
                        },
                data=payload)
r.json()

# Visualizar os identificadores
job_id ="3cf2faf653502b3b458d"

headers = {'Authorization': 'antoniosousa.js98@gmail.com bf6acfa27c682e8b136d'}
payload={'action':'logs'}
r = requests.get('https://www.cancergenomeinterpreter.org/api/v1/%s' % job_id, headers=headers, params=payload)
r.json()

# Download dos Resultados (file.zip)
job_id ="3cf2faf653502b3b458d"

headers = {'Authorization': 'antoniosousa.js98@gmail.com bf6acfa27c682e8b136d'}
payload={'action':'download'}
r = requests.get('https://www.cancergenomeinterpreter.org/api/v1/%s' % job_id, headers=headers, params=payload)
with open('/content/CGI/file.zip', 'wb') as fd:
    fd.write(r._content)

# Descompactar o arquivo "file.zip"
!unzip /content/CGI/file.zip -d /content/CGI/

# Resultado: "alterations.tsv"
import pandas as pd
pd.read_csv('/content/CGI/alterations.tsv', sep='\t', index_col=False, engine='python')

# Resultado: "alterations.tsv"
pd.read_csv('/content/CGI/biomarkers.tsv',sep='\t',index_col=False, engine= 'python')

# Deletar o Job do CGI
headers = {'Authorization': 'antoniosousa.js98@gmail.com bf6acfa27c682e8b136d'}
r = requests.delete('https://www.cancergenomeinterpreter.org/api/v1/%s' % job_id, headers=headers)
r.json()
```
