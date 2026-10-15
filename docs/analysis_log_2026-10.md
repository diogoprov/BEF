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
