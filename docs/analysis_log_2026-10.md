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

## 9. T99 (Downing 2005, Ecology) — passos 1–3 (`analysis/10_build_data_t99.R`, `11_models_t99.R`)

- 27 composições × 4 tanques; composição 0 = controle sem espécies (excluída).
  Células do desenho: composição de grupos funcionais × riqueza
  (M_6, G_6, MG_6, GP_6, MGP_6, MGP_18). 13 amostragens irregulares
  (dias 0–72); decomposição em 9.
- Mapeamento de códigos resolvido pelos dados; `Ecology2005.txt` tem 2 erros
  (comp. 5 sem *L. clamitans*; comp. 16 com *Buenoa* em vez de *Dineutus*),
  confirmados pela biomassa final. Fonte: planilha do desenho.
- Pulso ácido entre os períodos 4 (dia 19) e 5 (dia 22), inferido das médias
  [verificar em Downing & Leibold 2010]. Resistência = 5–6; recuperação = 7–13.
- Variância das métricas iniciais dentro das células: PD 1,8%, MPD 2,2%,
  VPD 3,8%, MNTD 16%.
- Heterocedasticidade temporal forte em prod/resp (DP varia 3–5×); primeira
  versão Student-t estimou nu ≈ 2. Versão final: gaussiana com sigma ~ ftime
  (prod, resp); Student-t (dc). Versão Student-t em `output/tables/sensitivity/`.
- (a) Dentro das células: nenhuma inclinação credível (98,33%).
- (b) Células vs métrica única: células explicam 89% (prod) e 79% (resp) da
  variância entre composições; VPD sozinho ~70%/68%, MPD 22–28%, PD ~0.
  VPD: prod 0,11 [0,06; 0,16], resp 0,12 [0,06; 0,18]. dc: nada.
- (c) Distúrbio: resp resistência ≈ −1,1 DP; recuperação 0,90 (S6), 1,13 (S18).
  MPD (inicial e final) modula negativamente a fase de recuperação de prod
  (−0,10; −0,08) e resp (−0,15; −0,13), todos fora do zero a 98,33%; mesmo
  sinal em T98 para resp (−0,10 [−0,26; 0,05]). Tentativo: checar influência
  de composições individuais.
- Convergência: R-hat ≤ 1,01 exceto 2 modelos com VPD (≤ 1,019) por trade-off
  entre efeitos fixos de célula e tempo (parâmetro focal ok). Corrigir com
  contrastes soma-zero em ftime no fechamento.

## 10. PD dentro dos grupos tróficos (`analysis/17_within_trophic_check.R`)

- Pergunta: a MPD calculada *dentro* de cada grupo trófico (macrófitas,
  herbívoros, predadores; composição estocada) varia entre composições dentro
  das células do desenho, e assim escapa do confundimento com a estrutura
  trófica? Apenas descritivo, sem modelos.
- Há variação, mas só com poucas espécies por grupo. CV da MPD dentro do
  grupo: 0,19–0,31 (T98, 3 spp.) e 0,19–0,60 (T99, 2–3 spp.), contra
  0,04–0,09 com 5–6 spp. (5 de 7 ou 6 de 8 espécies do pool: pouco espaço
  combinatório). O contraste útil se resume a ~7 composições por experimento
  (S9 no T98; MGP_6 no T99).
- MPD média dentro dos grupos: 72% da variância dentro dos níveis de riqueza
  no T98 (MPD da comunidade inteira: 11%), quase toda em S9 (DP 73 vs 12 Myr
  em S15); 9% dentro das células no T99 (DP 133 Myr em MGP_6).
- A variação vem de poucas bifurcações profundas, cada uma ligada a táxons
  reconhecíveis: girinos vs invertebrados (1565 Myr) e caramujos vs
  crustáceos/insetos (1172) nos herbívoros; Coleoptera vs Hemiptera (744) nos
  predadores; monocotiledôneas vs *Ceratophyllum*/*Myriophyllum*/*Utricularia*
  (~284–326) nas macrófitas. A MPD dentro de grupo funciona como índice de
  presença de *Lithobates*, caramujos e besouros → confundida com identidade
  de 1–2 espécies (mesmo problema do modelo do Box, em escala menor).
- Conclusão: a opção "PD dentro dos níveis tróficos" não dá um teste limpo.
  Usar como material suplementar para sustentar o enquadramento do artigo
  (PD multitrófica não se separa de estrutura trófica nem de identidade com
  estes desenhos); efeitos nulos dentro dos níveis tratados como não
  informativos. Título do abstract a revisar ("does not predict" excede o que
  o desenho sustenta).
- Pendência: as idades da árvore parecem altas (vertebrado–protostômio
  1565 Myr); conferir a calibração do pseudocronograma [não verificado].
- Próximo: Box de desenho experimental (simulação) — ver plano.

## 11. Box de desenho experimental (`analysis/18`–`20`)

- Pergunta: que desenho separaria um efeito de parentesco (θ, efeito da MPD
  dentro dos grupos) da estrutura trófica e da identidade de espécies, e
  quando σ²_phy (sinal filogenético nos efeitos de espécie) seria separável
  de σ²_id? Modelo gerador na parametrização do script 13 (kernels P C Pᵀ e
  P Pᵀ, projetados ortogonais aos efeitos fixos); escala de médias por
  composição; s_phy, s_id ∈ {0,15/0,15; 0,21/0; 0/0,21}, s_comp 0,10,
  s_tanque 0,30, 4 tanques; θ ∈ {0; 0,1; 0,2} DP por DP da MPD entre
  composições aleatórias do pool; 24/48/96 composições; 1000 réplicas
  (7,2 min, 17 núcleos).
- Desenhos: D1 montagem aleatória tipo Downing (k = 2 ou 6 por grupo); D2
  gradiente de MPD (pool T99, k = 3); D3 gradiente + incidência balanceada;
  D4b grupo focal = macrófitas T99 (8 spp.), resto fixo; D4 grupo focal com
  20 spp. em árvore pure-birth (10 árvores), resto fixo.
- Mecanismo: no D2 a MPD correlaciona 0,72 com a presença de
  *L. catesbeianus*; o balanceamento reduz para 0,41 (D3); D4 0,34.
- Poder (θ = 0,1, phy+id, modelo com kernels), 24 composições: D1 0,33; D2
  0,71; D3 0,86; D4 0,78; D4b 0,60. Com 96: todos ≥ 0,87.
- **Erro tipo I com composições tratadas como independentes (equivalente ao
  `(1 | comp)` iid) cresce com o número de composições**: 0,15–0,22 com 24 e
  0,26–0,52 com 96, porque composições compartilham espécies e o número
  efetivo de unidades satura no tamanho do pool. Com kernels: 0,10–0,13
  (24) e 0,06–0,08 (96); o excesso em 24 provavelmente vem da estimativa
  REML das variâncias na fronteira (checar com brms, script 19).
- **Implicação para os resultados do artigo**: os modelos com `(1 | comp)`
  (T98 21 e T99 26 composições) estão na faixa em que o erro tipo I nominal
  de 5% vira ~15–20%. Reforça a leitura dos poucos efeitos "credíveis"
  (ex.: MPD × recuperação no T99) como não confirmatórios.
- σ²_phy vs σ²_id: intervalo 10–90% da fração filogenética cobre ~0–1 em
  quase todos os desenhos e tamanhos; só D4 com 96 composições começa a
  separar (q10 0,53 quando a verdade é 1; q90 0,61 quando é 0). Separar
  efeitos de espécie filogeneticamente estruturados de identidade exige pool
  grande dentro do grupo e ~100 composições.
- Figura: `figures/box_design_simulation.{png,pdf}` (paleta checada para
  daltonismo; formas + rótulos diretos como codificação secundária).
- Ressalvas: θ definido por DP do pool favorece desenhos com seleção; D4 usa
  árvores pure-birth; D4b repete composições com 96 (só 56 triplets).
- Checagem bayesiana (script 19; D2 e D4, 24 composições, θ = 0, phy+id,
  100 réplicas cada; priors normal(0, 1) em b e Intercept, exponential(2) em
  sd e sigma; adapt_delta 0,99): IC 95% de θ exclui zero em 5% (D2) e 3% (D4)
  das réplicas (MCSE ~0,02), contra 12–13% do Wald com REML → o excesso do
  REML vem da estimativa pontual das variâncias.
- Diagnósticos por parâmetro (rodada repetida com as mesmas sementes): no D2
  os problemas estão em sd_id, sigma e sd_phy (R-hat máx. 1,19/1,17/1,14;
  ESS bulk mín. 17–24) — troca entre variância de identidade e residual com
  24 composições, coerente com a não identificabilidade do painel d. θ
  convergiu (R-hat < 1,01, ESS > 400) em 94/100 réplicas (D4: 99/100).
  Exclusão do zero: D2 5,3% (θ convergido, n = 94), 6,8% (tudo convergido,
  sem divergências, n = 59); D4 3,0% (n = 99) e 1,4% (n = 73). Todos
  compatíveis com 5% dado o MCSE (0,02–0,03). Conclusão: com a incerteza das
  variâncias propagada, o modelo com kernels mantém o erro tipo I nominal
  para θ; a inflação de 10–13% em 24 composições (painel c) é do Wald/REML.
  Prior mais apertada nos sd não testada: o problema está nas variâncias, que
  não são o alvo, e θ já está bem amostrado.
