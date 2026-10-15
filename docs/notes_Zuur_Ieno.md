# Notas da conversa com Alain Zuur e Elena Ieno (experimento Amy)

Transcrição das anotações originais:

- Tratar MPD como fator fixo (como fator), dar IDs diferentes; ver se funciona e
  se o modelo converge, e depois testar como contínuo.
- Usar o ID da composição como aleatório e o MPD como fixo.
- Primeiro testar se o tempo entra como aleatório; checar autocorrelação nos
  resíduos. Se não funcionar:
  - Time -> autoregressive AR1()?
  - `(1 + week | ID_pond)` => alternativa para tratar as repeated measures,
    random intercept
  - `AR1(0 + week | ID_pond)` => outra possibilidade mais complicada, último caso!
- P-value = 0.05/3 (Bonferroni) devido às múltiplas variáveis resposta =>
  alfa = 0.0167.
- Multivariate mixed modelling seria o correto, mas é demais; então tornar os
  p-values mais conservadores.

## Observação (Claude, 2026-10-15): identificabilidade

No `data_diogo.csv` há 21 composições e 21 valores únicos de PD (e MPD): cada
composição tem um valor único. Logo, **MPD como fator + `(1 | comp)` não é
identificável**: o fator MPD é a própria composição. Opções coerentes:

1. MPD contínuo + `(1 | comp)`: a variação entre composições não explicada pelo
   MPD vai para o efeito aleatório (é o teste correto do efeito do MPD, com
   21 unidades de replicação).
2. MPD como fator com `(1 | tank)` (sem `comp`), equivalente a testar
   composição.

O tanque (`tank`, 4 por composição) é a unidade das medidas repetidas e deve
entrar como `(1 | comp/tank)`. AR(1) no tempo dentro de tanque só se os
resíduos indicarem autocorrelação (em `brms`: `ar(time = time, gr = tank)`).
Alternativa ao Bonferroni: modelo multivariado em `brms`
(`mvbind(dc, resp, prod) ~ ...`), que é viável.
