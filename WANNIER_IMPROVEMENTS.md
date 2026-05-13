# Wannier Class Improvements & Sequential Pipeline Tasks

## Outstanding Issues

Certo! Aqui está o plano de refatoração que sugeri:

## 🎯 PLANO DE REFATORAÇÃO WANNIER_WORKFLOW.PY

### **PROBLEMA ATUAL:**
- Arquivo: 1341 linhas
- Método `run()`: 383 linhas (muito grande, difícil manter)
- Múltiplas responsabilidades misturadas

### **SOLUÇÃO PROPOSTA:**

#### **ETAPA 1 - Quebrar o método run() em sub-métodos (PRIORIDADE 1)**
Dividir `WannierWorkflow.run()` em 6 sub-métodos temáticos:
- [ ] `_run_scf_stage()` - executar SCF
- [ ] `_run_nscf_stage()` - executar NSCF  
- [ ] `_run_projwfc_stage()` - executar PROJWFC
- [ ] `_run_pw2wannier_stage()` - executar pw2wannier90
- [ ] `_run_wannier_stage()` - executar wannier90
- [ ] `_run_bands_stage()` - executar band structure (opcional)

Depois o `run()` fica apenas orquestrando esses 6 sub-métodos em sequência.

#### **ETAPA 2 - Extrair wannier_helpers.py (PRIORIDADE 2)**
Mover 8 funções standalone (459 linhas) para novo arquivo:
- [ ] `run_pw2wannier()`
- [ ] `run_projwfc()`
- [ ] `run_wannier90()`
- [ ] `generate_pw2wannier_input()`
- [ ] `generate_seedname_win()`
- [ ] `parse_projwfc_output()`
- [ ] `suggest_nbnd_from_pseudos()`
- [ ] `_make_queue_fallback()`

**Arquivo:** `xespresso/workflow/wannier_helpers.py`
**Depois atualizar imports** em wannier_workflow.py

#### **ETAPA 3 - Extrair wannier_analysis.py (PRIORIDADE 3)**
Mover análises e validação (145 linhas) para novo arquivo:
- [ ] `validate_wannier_quality()` (método)
- [ ] `compare_bands()` (método)
- [ ] Métodos auxiliares de validação

**Arquivo:** `xespresso/workflow/wannier_analysis.py`
**Depois atualizar imports** em wannier_workflow.py

#### **ETAPA 4 - Resultado Final**
- wannier_workflow.py ~ 400 linhas (core orquestração apenas)
- `wannier_helpers.py` ~ 459 linhas (helpers independentes)
- `wannier_analysis.py` ~ 145 linhas (análise/validação)
- **Total**: dividido em 3 arquivos focados, muito mais mantenível

### **ORDEM DE EXECUÇÃO:**
1. Etapa 1 (quebrar run() em sub-métodos) ← **COMEÇA AQUI**
2. Testar imports e funcionalidade
3. Etapa 2 (extrair helpers.py)
4. Testar novamente
5. Etapa 3 (extrair analysis.py)
6. Teste final completo
7. Commit: "refactor: modularize WannierWorkflow"

---

### 1. Wannier Class Feature Enhancements
**Location:** `xespresso/post/wannier90.py` / `xespresso/workflow/wannier_workflow.py`

**Missing Features/Parameters:**
- Additional Wannier90 input parameters (e.g., `use_bloch_phases`, `exclude_bands`, `bmin`, `bmax`, etc.)
- Better control over projection types and atomic projections
- Support for different guiding centers initialization methods
- Plotting and visualization parameters
- Restart/recovery from incomplete calculations
- Post-processing features (band structure interpolation, DOS calculation with Wannier)

**Action Items:**
- [ ] Expand `EspressoWannier90` class with more input parameters
- [ ] Add validation for parameter combinations
- [ ] Document all supported Wannier90 keywords with examples

### 2. Sequential Pipeline for SLURM Scheduler

**Context:** 
When users use SLURM scheduler (or similar) with multi-stage workflows (e.g., wannier90 -pp → pw2wannier90 → wannier90), the current implementation generates separate `srun` commands that may execute in parallel instead of sequentially.

**Problem:**
- Multiple `srun` calls in a single job script can be interpreted as independent parallel processes
- No guaranteed sequential execution without explicit dependencies
- Current solutions tested and rejected:
  - `set -e` flag: Not acceptable (exit code control not desired)
  - `bash -c` wrapper with single srun: Not acceptable (preference for simpler format)

**Current State:**
- Implementation: Multiple separate `srun` calls without guarantees on execution order
- This is a **blocking issue** for production use with SLURM scheduler and multi-step workflows

**Potential Solutions to Explore:**
1. **Job Dependencies (SLURM):** Generate 3 separate job scripts with `#SBATCH --dependency=afterok:job_id` directives
   - Pros: Clean, SLURM-native, explicit control
   - Cons: Requires manual job tracking or wrapper script
   
2. **Custom Scheduler Wrapper:** Create an abstract scheduler pipeline that handles multi-step workflows
   - Pros: Scheduler-agnostic, integrates with xespresso architecture
   - Cons: Significant refactoring needed
   
3. **Workflow Engine Integration:** Use a proper workflow engine (Snakemake, Nextflow)
   - Pros: Production-ready, handles complex pipelines
   - Cons: External dependency, different execution model
   
4. **Documentation & User Guidelines:** Document that users should either:
   - Run stages individually with `sbatch` and wait for completion
   - Use job submission script that monitors completion before submitting next stage
   - Pros: No code changes, keeps simplicity
   - Cons: Manual process, error-prone

**Affected Code:**
- `xespresso/workflow/wannier_workflow.py` - `_generate_wannier_job_script()` method
- `xespresso/scheduler/` - scheduler abstraction layer
- `xespresso/post/` - post-processing tools (pw2wannier90, wannier90)

**Decision Needed:**
Determine which solution best fits the xespresso philosophy and project maturity level.

## Session Context

**Date:** May 13, 2026

**What was accomplished:**
1. Verified WannierWorkflow generates proper SLURM headers for all job_files (SCF, projwfc, NSCF, wannier_up, wannier_dn)
2. Identified that machine configuration (scheduler) properly propagates from CalculationWorkflow to post-processing tools
3. Discovered issue with sequential execution of 3-stage Wannier workflow when using SLURM scheduler
4. Explored multiple solutions (set -e, bash wrapper) but user rejected in favor of simpler format
5. Created this document to track improvements for future implementation

**Ready for Next Session:**
- WannierWorkflow SLURM integration is complete and tested
- Job_files are properly formatted for SLURM submission
- Sequential pipeline issue documented and awaiting design decision
