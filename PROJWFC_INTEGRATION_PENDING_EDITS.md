# STATUS: WANNIER90 WORKFLOW PROJWFC INTEGRATION - PENDING EDITS

## Resumo do Que Foi Feito

✓ Função `run_projwfc()` criada (linha 118)
✓ Função `parse_projwfc_output()` criada (linha 212)
✓ Exemplos e documentação criados
✓ Testes criados

## ❌ O QUE AINDA FALTA (Edições Pendentes)

### 1. Reordenar o Pipeline do Método `run()` 

**Ordem atual ERRADA:**
```
SCF → Bands → NSCF → pw2wannier90 → wannier90
```

**Ordem desejada CORRETA:**
```
SCF → Bands → PROJWFC → NSCF → pw2wannier90 → wannier90
```

---

### EDIT 1: Adicionar PROJWFC entre Bands e NSCF

**Localização:** Após a seção "# ====== STAGE 2 (OPTIONAL): BAND STRUCTURE ======" (por volta da linha 550-570)

**O que fazer:**
Inserir a seguinte seção ANTES da seção "# ====== STAGE 3: NSCF ======"

```python
        # ====== STAGE 3 (OPTIONAL): PROJWFC ======
        projwfc_projections_suggested = None
        if run_projwfc_analysis:
            print(f"\n[{stage_count}/{total_stages}] Running projwfc (projection analysis for Wannier guidance)...")
            stage_count += 1
            try:
                # Use SCF results for projwfc
                run_dir = str(Path(scf_calc.directory).resolve())
                prefix = scf_calc.prefix
                
                projwfc_result = run_projwfc(
                    run_dir=run_dir,
                    prefix=prefix,
                    blocking=blocking,
                    queue=self.queue
                )
                self.results['projwfc'] = projwfc_result
                
                if projwfc_result['status'] in ['finished', 'finished_with_warnings', 'submitted']:
                    print(f"✓ projwfc completed: {projwfc_result['status']}")
                    if 'outputs' in projwfc_result and projwfc_result['outputs']['pdos']:
                        print(f"  PDOS analysis: {Path(projwfc_result['outputs']['pdos']).name}")
                        
                        # Try to extract suggestions from PDOS output
                        suggestions = parse_projwfc_output(run_dir, prefix)
                        if suggestions:
                            self.results['projwfc_analysis'] = suggestions
                            projwfc_projections_suggested = suggestions.get('projections', None)
                            if projwfc_projections_suggested:
                                print(f"  Suggested projections from PDOS: {projwfc_projections_suggested}")
                        
                        print(f"  ➜ Check {run_dir}/{prefix}.pdos* for detailed orbital contributions")
                else:
                    print(f"✗ projwfc {projwfc_result['status']}: {projwfc_result.get('message', 'Unknown error')}")
                    print(f"  Continuing with user-specified projections...")
                    self.results['projwfc'] = None
            except Exception as e:
                print(f"⚠ projwfc failed (non-critical): {e}")
                print(f"  Continuing with user-specified projections...")
                self.results['projwfc'] = None
```

---

### EDIT 2: Renumerar os estágios NSCF, pw2wannier90 e wannier90

Precisa mudar:
- "STAGE 3: NSCF" → "STAGE 4: NSCF"
- "STAGE 4: pw2wannier90" → "STAGE 5: pw2wannier90"  
- "STAGE 5: wannier90" → "STAGE 6: wannier90"

---

### EDIT 3: Atualizar contagem de estágios

**Localização:** Por volta da linha 525-530

**Atual:**
```python
        total_stages = 6 if (run_bands_validation and run_projwfc_analysis) else (5 if run_bands_validation else 4)
```

**Mudar para:**
```python
        total_stages = 5  # Base: SCF + NSCF + pw2wannier + wannier90
        if run_bands_validation:
            total_stages += 1  # Add bands
        if run_projwfc_analysis:
            total_stages += 1  # Add projwfc
```

---

### EDIT 4: Adicionar método `get_projwfc_analysis()` à classe

**Localização:** Após o método `get_band_structure_calculator()` (por volta da linha 750)

**O que adicionar:**
```python
    def get_projwfc_analysis(self) -> Optional[Dict]:
        """Get PROJWFC analysis results.
        
        Returns:
            Dict with PDOS file paths and analysis, or None if PROJWFC was not run
        """
        if not self.results.get('projwfc'):
            return None
        
        projwfc_result = self.results['projwfc']
        analysis = {
            'status': projwfc_result.get('status'),
            'pdos_file': projwfc_result.get('outputs', {}).get('pdos'),
            'run_dir': projwfc_result.get('run_dir'),
            'job_id': projwfc_result.get('job_id'),
        }
        
        if 'projwfc_analysis' in self.results:
            analysis['parsed'] = self.results['projwfc_analysis']
        
        return analysis
```

---

### EDIT 5: Atualizar validate_wannier_quality() 

**Localização:** Método `validate_wannier_quality()` (por volta da linha 755+)

**Adicionar após 'has_band_structure':**
```python
            'has_projwfc': self.results.get('projwfc_available', False),
```

**Adicionar após as recomendações de band_structure:**
```python
        if validation['has_projwfc']:
            validation['recommendations'].append(
                "PROJWFC analysis available - Review PDOS to understand orbital contributions"
            )
        else:
            validation['recommendations'].append(
                "Consider re-running with run_projwfc_analysis=True for orbital analysis guidance"
            )
```

---

### EDIT 6: Atualizar seção de conclusão

**Localização:** Por volta da linha 720+, na seção de conclusão do método `run()`

**Atualizar a linha:**
```python
        self.results['band_structure_available'] = run_bands_validation and self.results['bands'] is not None
```

**Para incluir:**
```python
        self.results['band_structure_available'] = run_bands_validation and self.results['bands'] is not None
        self.results['projwfc_available'] = run_projwfc_analysis and self.results.get('projwfc') is not None
```

---

### EDIT 7: Atualizar mensagem de conclusão

**Localização:** Seção de mensagens finais do método `run()` (por volta da linha 730+)

**Adicionar após seção de Band Structure Validation:**
```python
        if run_projwfc_analysis and self.results.get('projwfc') is not None:
            print(f"\nProjection Analysis (PROJWFC):")
            print(f"  - PDOS files: {prefix}.pdos*")
            print(f"  - Shows orbital contributions to electronic structure")
            print(f"  - Use to refine Wannier projections in next iterations")
```

---

## Resumo das Edições Necessárias

| # | Tipo | Localização | Prioridade |
|---|------|------------|-----------|
| 1 | INSERT | Após Bands, antes NSCF | ALTA |
| 2 | REPLACE | Números dos stages | ALTA |
| 3 | REPLACE | Cálculo total_stages | MÉDIA |
| 4 | INSERT | Novo método get_projwfc_analysis() | ALTA |
| 5 | UPDATE | Método validate_wannier_quality() | MÉDIA |
| 6 | UPDATE | Armazenar projwfc_available | MÉDIA |
| 7 | UPDATE | Mensagens finais | BAIXA |

---

## Arquivo Já Criados (Completos)

✓ `/home/vinicius/scratch/projects/spresso/xespresso/workflow/wannier_workflow.py` - com funções run_projwfc() e parse_projwfc_output()
✓ `/home/vinicius/scratch/projects/spresso/examples/wannier_workflow_with_projwfc_example.py`
✓ `/home/vinicius/scratch/projects/spresso/docs/PROJWFC_IN_WANNIER_WORKFLOW.md`
✓ `/home/vinicius/scratch/projects/spresso/docs/PROJWFC_INTEGRATION_SUMMARY.md`
✓ `/home/vinicius/scratch/projects/spresso/tests/test_projwfc_integration.py`

---

## Próximos Passos

1. Fazer as 7 edições acima no arquivo wannier_workflow.py
2. Executar: `conda run -n spresso python tests/test_projwfc_integration.py`
3. Todos os 6 testes devem passar (atualmente 5/6 passam)
4. Pronto! Integração completa.

---

**Data:** 12 de Maio de 2026
**Status:** ⚠️ AGUARDANDO EDIÇÕES MANUAIS
