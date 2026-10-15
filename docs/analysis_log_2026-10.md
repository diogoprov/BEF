# Registro da análise — outubro de 2026

## 1. Recálculo do PD (`analysis/00_build_data.R`)

- Fonte bruta: aba "ecosystem functioning T98 long" da planilha da Amy.
  `data_diogo.csv` é a mesma tabela arredondada, com um erro (`dist = 6`
  no tanque 5, período 6).
- `treat.txt` = composição final; validado: é exatamente a união das espécies
  com biomassa final > 0 nos 4 tanques de cada composição (aba "species
  biomass T98").
- PD recalculado com a árvore re-datada (`phylog_dated_2002.txt`). Diferença
  de +21,1 Myr em 9 composições em relação ao PD antigo (árvore anterior +
  `drop.tip`). PSV da composição 6 corrigida (0,620 → 0,762).
- **Bug no código antigo**: `rep(mpd, each = 12)` com 24 linhas por composição
  → MPD/MNTD errados em 480 de 504 linhas em todas as análises de MPD/MNTD de
  2014–2019 (Tables 2, 4, 6).
- Sensibilidade aos comprimentos de ramo: r(bladj, Grafen) = 0,99 (PD),
  0,99 (MPD), 0,97 (MNTD), 0,99 (PSV).
- Correlações entre preditores (21 composições): PD–div 0,88; MNTD–div −0,72;
  MPD–div 0,14; MPD–PSV 1,00 (PSV é redundante).
- Métricas por tanque (`metrics_tank.csv`) estão disponíveis para análise de
  sensibilidade; 2 tanques têm só 1 espécie (MPD/MNTD indefinidos).

## 2. Modelos (`analysis/01_models.R`)

`y ~ 0 + div + div:zX + ftime + dist + dist:post + (1|comp) + (1|tank) + ar(time, gr = tank)`

- Respostas padronizadas: `prod` (escala original), `log(resp)`, `dc`
  (Student-t). Escolhas guiadas por PPC: log(prod) era assimétrico à esquerda;
  dc tem caudas pesadas.
- Priors: b ~ N(0, 1); sd ~ Exp(2); sigma ~ Exp(1); ar ~ N(0, 0,5);
  nu ~ Gamma(2, 0,1); cor/rescor ~ LKJ(2).
- AR(1): ar[1] = 0,30 (prod), 0,52 (resp), 0,05 (dc). LOO favorece AR(1) para
  resp (Δelpd = −20,2 ± 9,9) e prod (−4,9 ± 5,9); dc indiferente. Mantido em
  todos os modelos.
- Diagnósticos (`output/tables/diagnostics.csv`): 0 divergências; R-hat ≤ 1,01
  em todos os parâmetros focais. Exceções em parâmetros secundários:
  `cor_tank(prod, dc)` no multivariado (R-hat 1,017; média −0,01) e `lp__` em
  dc ~ MNTD (1,013; 4 transições no treedepth máximo).

## 3. Resultados principais (escala padronizada, DP da resposta por DP do preditor)

- **Nenhum efeito credível de PD, MPD ou MNTD** sobre produtividade,
  respiração ou decomposição; todos os IC 95% (e 98,33%) incluem zero.
  Ex.: inclinação média do MPD — prod 0,12 [−0,19; 0,44], resp 0,12
  [−0,23; 0,45], dc 0,18 [−0,23; 0,60]. O modelo multivariado concorda.
- Inclinações em S15 são muito incertas: o MPD varia pouco entre as
  composições de 15 espécies.
- Pulso ácido (dist × post): prod −0,45 [−0,68; −0,22]; resp −0,92
  [−1,12; −0,71]; dc −0,29 [−0,47; −0,11].

## 4. Pendências

- Os resultados antigos de efeito quadrático do PD e de MPD vieram de `lm`
  tratando as 504 observações como independentes (pseudorreplicação) e/ou do
  MPD embaralhado; não devem ser reportados.
- Discutir: PD a partir da composição final (resposta parcial ao experimento)
  e não da inicial manipulada.
- Opcional: modelos com métricas por tanque; interação MPD × distúrbio
  (estabilidade).

## 5. Composição inicial (suplemento da *Nature*)

- `analysis/00a_initial_composition.py` extrai a Table 1 do SI de Downing &
  Leibold (2002) (`refs/`, não versionado): "–" não estocada; "0" estocada e
  extinta; valor = estocada e sobreviveu; itálico = invasora. Validação: todas
  as 21 composições têm 1/3/5 espécies por grupo funcional; sobreviventes +
  invasoras reproduzem exatamente o `treat.txt`. Invasoras: *Neoplea* (7
  composições), *Physa* (7), *Trichocorixa* (3), *Belostoma* (1).
- Métricas da composição inicial usam a árvore de 21 espécies (commit
  `4b4da4d`, `data/raw/phylogeny/phylog_dated_2002_21spp_4b4da4d.txt`).
- **As métricas iniciais são quase determinadas pela riqueza**: fração da
  variância dentro dos níveis de riqueza = 2% (PD), 11% (MPD), 86% (VPD).
  Causa: cada composição tem proporções iguais de macrófitas, herbívoros e
  predadores, e a divergência planta–animal domina as distâncias.

## 6. Passo 2 (`analysis/02_initial_composition.R`)

- PD, MPD, VPD fixados a priori (Tucker et al. 2017), padronizados dentro de
  cada nível de riqueza; composição inicial e final.
- Nenhuma inclinação credível. Maior: dc ~ MPD inicial 0,25 [−0,06; 0,57].
- Variância entre composições pequena e mal estimada (sd_comp: prod 0,16,
  resp 0,20, dc 0,29; sigma 0,37–0,58); R² de composição não informativo.
- Riqueza: S15 ≈ +0,4 DP acima de S9 em prod e resp; S3 ≈ S9.

## 7. Passo 3 (`analysis/03_disturbance.R`)

- Resistência (período 3) e fase de recuperação (4–6) por nível de riqueza.
  Respiração: queda ≈ −1 DP em todos os níveis; recuperação (late − early)
  S3 0,42 [0,12; 0,74], S9 0,18 [−0,13; 0,49], S15 0,71 [0,39; 1,02].
- Modulação filogenética: nenhuma credível em 98,33%. Único IC 95% sem zero:
  resp, PD inicial × fase de recuperação −0,19 [−0,35; −0,04] (1 de 24
  testes; esperado ~1 por acaso).

## 8. Passo 4 — avaliação

- Modelo de identidade de espécie como efeito aleatório de composição com
  covariância K = P C Pᵀ (filogenético) e P Pᵀ (identidade), kernels
  duplo-centrados; ligação com DI models (Kirwan et al. 2009; Connolly et al.
  2011): com proporções iguais, a interação dependente da distância
  filogenética se reduz a θ(S−1)/(2S)·MPD.
- T98: alinhamento entre kernels 0,87; posto 18; variância entre composições
  pequena → separação filogenia vs identidade fracamente identificada.
- T99 (Downing 2005, Ecology): 26 composições com composição de grupos
  funcionais variável; alinhamento 0,70; razão de traços Kphy/Kid 1,57
  (T98: 0,26). Candidato melhor para o passo 4 e para replicação.
