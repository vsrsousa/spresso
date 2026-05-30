# Plano de Trabalho: Implementação do Módulo EOS (Equação de Estado)

**Data de Criação**: 29 de Maio de 2026  
**Objetivo**: Implementar funcionalidade completa de Equação de Estado (Birch-Murnaghan) no xespresso  
**Escopo**: Novo módulo `xespresso/workflow/eos_workflow.py` com integração ao sistema existente

---

## 📋 Estrutura do Plano

### Fase 1: Análise e Design
**Status**: ⏳ Não iniciado  
**Objetivo**: Definir arquitetura detalhada

#### 1.1 Análise de Requisitos
- [ ] Revisar implementações existentes em ASE (EquationOfState)
- [ ] Verificar dependências necessárias (scipy, numpy, matplotlib)
- [ ] Mapear integração com CalculationWorkflow
- [ ] Definir interface pública (métodos, parâmetros)
- **Saídas**: Documento de requisitos, definição de API

#### 1.2 Design Arquitetural
- [ ] Definir classe EOSWorkflow e herança (Base ou CalculationWorkflow?)
- [ ] Desenhar fluxo de dados E-V
- [ ] Especificar estrutura de resultados
- [ ] Planejar paralelização de SCFs
- **Saídas**: Diagrama de arquitetura, interface de classes

#### 1.3 Planejamento de Testes
- [ ] Definir casos de teste (Si, Fe, Au)
- [ ] Planejar dados mock para testes rápidos
- [ ] Especificar critérios de validação
- **Saídas**: Estratégia de testes

---

### Fase 2: Criação do Arquivo Base
**Status**: ⏳ Não iniciado  
**Objetivo**: Criar estrutura inicial com todas as classes

#### 2.1 Criar `xespresso/workflow/eos_workflow.py`
- [ ] Adicionar imports (logging, numpy, scipy, matplotlib)
- [ ] Criar classe EOSWorkflow (estrutura básica)
- [ ] Adicionar `__init__` com parâmetros principais
- [ ] Implementar logging
- [ ] Adicionar docstrings completas
- **Artefato**: `xespresso/workflow/eos_workflow.py` (versão 1)

#### 2.2 Definir Constantes e Configurações
- [ ] Volumes padrão (range: 95%-105%)
- [ ] Número de pontos (padrão: 7)
- [ ] Tolerâncias para fitting
- [ ] Configurações de paralelização
- **Artefato**: Seção CONSTANTS em eos_workflow.py

---

### Fase 3: Implementar Volume Scaling
**Status**: ⏳ Não iniciado  
**Objetivo**: Ferramentas para escalar volume da estrutura

#### 3.1 Função `scale_atoms_isotropic(atoms, scale_factor)`
- [ ] Escalona volume por fator (1.0 = original)
- [ ] Mantém ângulos da célula
- [ ] Ajusta posições internas proporcionalmente
- [ ] Testes: verificar volume e distâncias relativas
- **Saída**: Método scale_volume_uniformly()

#### 3.2 Função `generate_volume_range(atoms, volume_range, n_points)`
- [ ] Gera lista de volumes uniformemente espaçados
- [ ] Volume_range: (min_factor, max_factor) ex: (0.95, 1.05)
- [ ] Retorna lista de factors para scaling
- [ ] Testes: verificar distribuição uniforme
- **Saída**: Método generate_volume_range()

#### 3.3 Função `create_eos_structures(atoms, factors)`
- [ ] Cria estruturas escaladas para cada fator
- [ ] Retorna dicionário {factor: atoms_scaled}
- [ ] Valida estruturas (distâncias mínimas, etc)
- **Saída**: Método create_scaled_structures()

---

### Fase 4: Implementar Coleta de Dados E-V
**Status**: ⏳ Não iniciado  
**Objetivo**: Executar SCFs e coletar resultados

#### 4.1 Função `run_scf_for_volume(atoms, scale_factor, label, protocol)`
- [ ] Reutiliza CalculationWorkflow.run_scf()
- [ ] Cria diretório com label contendo factor
- [ ] Executa e retorna (E, V)
- [ ] Captura erros graciosamente
- **Saída**: Método _run_single_eos_point()

#### 4.2 Função `run_eos_study(volume_range, n_points, protocol, parallel)`
- [ ] Gera estruturas escaladas
- [ ] Executa SCFs (paralelo se True)
- [ ] Coleta E-V pairs
- [ ] Ordena por volume
- [ ] Retorna pandas DataFrame com colunas: [volume, energy, factor]
- **Saída**: Método run_eos_study()

#### 4.3 Integração com Paralelização
- [ ] Use ThreadPoolExecutor ou ProcessPoolExecutor
- [ ] Controle máximo de workers baseado em n_points
- [ ] Implementar progress tracking
- [ ] Tratamento de exceções robustos
- **Saída**: Método _run_parallel_eos_points()

---

### Fase 5: Implementar Birch-Murnaghan Fitting
**Status**: ⏳ Não iniciado  
**Objetivo**: Ajustar dados aos modelos EOS

#### 5.1 Função `birch_murnaghan_fit(volumes, energies)`
- [ ] Implementar EOS de Murnaghan de 2ª ordem
- [ ] E(V) = E_0 + B_0*V_0/(B_0'-1) * [(V_0/V)^(B_0'-1)/(B_0'-1) + 1]
- [ ] Usar scipy.optimize.minimize ou curve_fit
- [ ] Retornar {E_0, V_0, B_0, B_0_prime}
- [ ] Validar convergência
- **Saída**: Método fit_birch_murnaghan()

#### 5.2 Função `fit_murnaghan_alt(volumes, energies)` (alternativa)
- [ ] Versão simplificada (quadrática ou polinômio 3º grau)
- [ ] Para comparação com BM
- [ ] Útil se BM não convergir
- **Saída**: Método fit_murnaghan_simplified()

#### 5.3 Validação do Fit
- [ ] Calcular R² (coeficiente determinação)
- [ ] Verificar resíduos
- [ ] Alertar se R² < 0.99
- [ ] Retornar diagnóstico
- **Saída**: Método validate_fit()

---

### Fase 6: Implementar Análise e Extração de Propriedades
**Status**: ⏳ Não iniciado  
**Objetivo**: Calcular e retornar propriedades físicas

#### 6.1 Função `get_eos_properties()`
- [ ] Retornar dicionário com:
  - V_0: volume de equilíbrio (Ų)
  - E_0: energia em V_0 (eV)
  - B_0: bulk modulus (GPa)
  - B_0_prime: primeira derivada do bulk modulus
  - Energy_shift: (E_0 - E_min_SCF) para referência
- [ ] Incluir erros/uncertainties se possível
- **Saída**: Método get_eos_properties()

#### 6.2 Função `predict_energy_at_volume(V_new)`
- [ ] Use parâmetros BM para prever E(V_novo)
- [ ] Útil para interpolar/extrapolar
- **Saída**: Método predict_energy()

#### 6.3 Função `calculate_pressure(volumes, energies, V_new)`
- [ ] Calcular pressão: P = -dE/dV (usando BM)
- [ ] Retornar em GPa
- **Saída**: Método calculate_pressure()

---

### Fase 7: Implementar Visualização (Plotting)
**Status**: ⏳ Não iniciado  
**Objetivo**: Gerar gráficos de E-V com fit

#### 7.1 Função `plot_eos_curve(save_path=None)`
- [ ] Plot dados E-V como scatter
- [ ] Plot curva BM ajustada como linha
- [ ] Marcar V_0 com linha vertical
- [ ] Legenda, eixos com unidades
- [ ] Salvar em PNG se save_path fornecido
- **Saída**: Método plot_eos_curve()

#### 7.2 Função `plot_residuals(save_path=None)`
- [ ] Plot (E_data - E_BM) vs V
- [ ] Mostrar distribuição de erros
- [ ] Útil para diagnóstico
- **Saída**: Método plot_residuals()

#### 7.3 Função `plot_pressure(save_path=None)`
- [ ] Plot P(V) derivado de BM
- [ ] Mostrar P=0 (equilíbrio)
- **Saída**: Método plot_pressure()

---

### Fase 8: Criar Testes Unitários
**Status**: ⏳ Não iniciado  
**Objetivo**: Testes robustos de todos componentes

#### 8.1 Criar `tests/test_eos_workflow.py`
- [ ] Teste estrutura: test_*.py com fixtures
- [ ] Mock calculators quando necessário
- [ ] Usar pytest como framework

#### 8.2 Testes de Volume Scaling
- [ ] test_scale_atoms_isotropic()
- [ ] test_generate_volume_range()
- [ ] test_create_eos_structures()
- [ ] Verificar: volumes corretos, estrutura preservada, distâncias mantidas

#### 8.3 Testes de Fitting
- [ ] test_birch_murnaghan_fit() com dados sintéticos
- [ ] test_fit_quality_high_noise()
- [ ] test_fit_warnings()
- [ ] Verificar R² > 0.99 para dados bons

#### 8.4 Testes de E2E
- [ ] test_run_eos_study_mock()
- [ ] Mock CalculationWorkflow.run_scf()
- [ ] Verificar coleta e ordem de dados
- [ ] Testes com paralelização

#### 8.5 Integração com Fixtures ASE
- [ ] Fixture: bulk_si, bulk_fe, bulk_au
- [ ] Fixture: pseudopotentials_dict
- [ ] Fixture: protocol_presets

---

### Fase 9: Integração no Módulo Principal
**Status**: ⏳ Não iniciado  
**Objetivo**: Exportar e disponibilizar para usuários

#### 9.1 Atualizar `xespresso/workflow/__init__.py`
- [ ] Adicionar import: `from xespresso.workflow.eos_workflow import EOSWorkflow`
- [ ] Adicionar a `__all__`: `'EOSWorkflow'`
- [ ] Testar import: `from xespresso import EOSWorkflow`

#### 9.2 Atualizar `xespresso/__init__.py`
- [ ] Adicionar EOS ao _EXPORTS se necessário
- [ ] Garantir que EOSWorkflow está acessível top-level

#### 9.3 Atualizar README.md
- [ ] Adicionar seção "EOS Workflow" com exemplo simples
- [ ] Mencionar em Features

---

### Fase 10: Criar Exemplo de Uso
**Status**: ⏳ Não iniciado  
**Objetivo**: Documentação e exemplos práticos

#### 10.1 Criar `examples/eos_workflow_example.py`
- [ ] Exemplo 1: EOS simples com Si
- [ ] Exemplo 2: EOS com Fe (ferro)
- [ ] Exemplo 3: Comparar protocolos (fast vs accurate)
- [ ] Exemplo 4: Carregar e reanalisar resultados salvos
- [ ] Comentários explicativos em português

#### 10.2 Criar documentação
- [ ] Docstring completa da classe
- [ ] Docstring para cada método público
- [ ] Notas sobre performance e paralelização
- [ ] Referências bibliográficas (Birch-Murnaghan)

---

### Fase 11: Verificação Final e Testes
**Status**: ⏳ Não iniciado  
**Objetivo**: Validar tudo funciona corretamente

#### 11.1 Testes de Aceitação
- [ ] Executar pytest em todos testes
- [ ] Cobertura mínima: 85%
- [ ] Sem warnings ou deprecations

#### 11.2 Testes de Integração
- [ ] Executar exemplo completo com dados reais
- [ ] Verificar saída em arquivo/console
- [ ] Validar plots gerados

#### 11.3 Revisão de Código
- [ ] PEP8 compliance
- [ ] Docstrings presentes
- [ ] Erros tratados graciosamente
- [ ] Logging apropriado

#### 11.4 Documentação Final
- [ ] README atualizado
- [ ] Docstrings completos
- [ ] Exemplos funcionando

---

## 📊 Rastreamento de Progresso

| Fase | Tarefa | Status | Responsável | Início | Fim |
|------|--------|--------|------------|--------|-----|
| 1 | Análise e Design | ⏳ | - | - | - |
| 2 | Arquivo Base | ⏳ | - | - | - |
| 3 | Volume Scaling | ⏳ | - | - | - |
| 4 | Coleta E-V | ⏳ | - | - | - |
| 5 | Fitting BM | ⏳ | - | - | - |
| 6 | Análise | ⏳ | - | - | - |
| 7 | Plotting | ⏳ | - | - | - |
| 8 | Testes | ⏳ | - | - | - |
| 9 | Integração | ⏳ | - | - | - |
| 10 | Exemplos | ⏳ | - | - | - |
| 11 | Verificação | ⏳ | - | - | - |

---

## 🎯 Critérios de Conclusão

### Por Fase
1. ✅ Design aprovado
2. ✅ Arquivo compilável sem erros
3. ✅ Todos métodos de scaling testados
4. ✅ E-V curve coletada com sucesso
5. ✅ BM fit com R² > 0.99
6. ✅ Propriedades calculadas corretamente
7. ✅ Plots geram sem erros
8. ✅ Cobertura de testes > 85%
9. ✅ Imports funcionam
10. ✅ Exemplos executam com sucesso
11. ✅ Tudo funciona ponta-a-ponta

### Global
- ✅ Nenhum erro em pytest
- ✅ Documentação completa
- ✅ Exemplo rodando
- ✅ Integrado em `xespresso/__init__.py`
- ✅ Backward compatible

---

## 📚 Referências

### Equação de Estado Birch-Murnaghan
- Birch, F. (1947): "Finite strain isotropic analysis"
- Fei, Y. (1995): "Thermal expansion", Mineral Physics & Crystallography

### Implementações Existentes
- ASE: `ase.eos.EquationOfState`
- pymatgen: `pymatgen.analysis.eos`

---

## 🔄 Próximas Etapas

1. **Após aprovação deste plano**: Iniciar Fase 1 (Análise)
2. **Estimativa de tempo**: 6-8 horas para implementação completa
3. **Paralelização possível**: Fases 3, 4, 7 podem ser feitas em paralelo após Fase 2

---

**Status Geral**: 📋 Planejado e Pronto para Iniciar
