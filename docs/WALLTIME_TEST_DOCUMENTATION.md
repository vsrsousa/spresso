# Walltime Job File Generation Test

## Overview

O teste `test_walltime_jobfile_mock.py` demonstra que o parâmetro `walltime` flui corretamente através de toda a cadeia workflow até ser escrito no arquivo jobfile SLURM.

## Fluxo Completo

```
run_slab_convergence(walltime='2:00:00')
         ↓
queue = {'resources': {'time': '2:00:00'}}
         ↓
Espresso.write_input()
         ↓
SlurmScheduler.write_script()
         ↓
Job File: #SBATCH --time=2:00:00
```

## Resultados do Teste

### ✅ Test 1: Walltime Reaches Job File
```
🔧 Queue Configuration:
   - scheduler: slurm
   - resources['time']: 2:00:00

📋 Job File Content:
#!/bin/bash

#SBATCH --job-name=test_walltime
#SBATCH --output=test_walltime.out
#SBATCH --error=test_walltime.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=2:00:00          ← ✓ Walltime presente!
#SBATCH --partition=gpu

pw.x    -in  test_walltime.pwi  >  test_walltime.pwo
```

### ✅ Test 2: Multiple Walltime Formats
Validado que diferentes formatos funcionam:
- `1:00:00` (1 hora)
- `2:30:45` (2.5 horas)
- `12:00:00` (12 horas)
- `00:30:00` (30 minutos)

Todos escrevem corretamente como `#SBATCH --time=<valor>` no jobfile.

### ✅ Test 3: Walltime vs Job_Timeout (Parâmetros Independentes)
```
Especificado:
   - walltime='2:00:00'      → Controla tempo máximo do job SLURM
   - job_timeout=7200        → Controla tempo de espera do Python

Resultado:
   ✓ Walltime no job file: #SBATCH --time=2:00:00
   ✓ Job_timeout armazenado: queue['job_timeout'] = 7200
   ✓ Ambos funcionam independentemente como projetado
```

### ✅ Test 4: Complete Workflow Integration
Demonstra a integração completa:
1. Parâmetros entram na configuração de queue
2. SlurmScheduler lê `queue['resources']`
3. Escreve `#SBATCH --time=...` no jobfile
4. Todos os parâmetros SLURM aparecem corretamente

## Como Usar o Teste

### Executar todos os testes:
```bash
cd /home/vinicius/projects/spresso
python tests/test_walltime_jobfile_mock.py
```

### Executar teste específico:
```python
from tests.test_walltime_jobfile_mock import test_walltime_reaches_jobfile
test_walltime_reaches_jobfile()
```

## O que o Mock Faz

O teste usa `@patch('xespresso.scheduler.check_slurm_available')` para:

1. **Simular SLURM instalado**: Sem precisar ter SLURM realmente instalado
2. **Testar em qualquer máquina**: CI/CD, laptops, etc
3. **Validar jobfile generation**: Verifica que o arquivo é escrito corretamente
4. **Rápido**: Não espera por jobs reais

## Exemplo de Uso Real

Quando você chamar no seu código:

```python
from xespresso.workflow.slab_workflow import SlabWorkflow

slab_wf = SlabWorkflow(...)

# Convergência com walltime e job_timeout
results = slab_wf.run_slab_convergence(
    surface_index=(1, 1, 1),
    vacuum_test=[5, 7, 9, 11],
    nlayers_test=[3, 4, 5],
    machine='medusa',
    code_version='7.4.1',
    
    # ← Parâmetros novos:
    walltime='2:00:00',      # Tempo máximo do job SLURM
    job_timeout=7200,        # Python espera máximo 2 horas
)
```

Internamente, isso:
1. Cria `queue['resources']['time'] = '2:00:00'`
2. Cria `queue['job_timeout'] = 7200`
3. Passa para `CalculationWorkflow(..., queue=queue)`
4. Espresso gera jobfile com `#SBATCH --time=2:00:00`
5. Job é submetido a SLURM com tempo máximo de 2 horas

## Validações Realizadas

O teste verifica:

✅ **Bash shebang**: `#!/bin/bash` presente
✅ **SBATCH directives**: `#SBATCH --...` presentes
✅ **Walltime directive**: `#SBATCH --time=<valor>` corretamente escrito
✅ **Outros parâmetros**: nodes, ntasks-per-node, partition, etc
✅ **Comando de execução**: `pw.x -in ... > ...` presente
✅ **Múltiplos formatos**: Diferentes valores de walltime validados
✅ **Integração completa**: Parâmetros fluem através de toda cadeia

## Scheduler Compatibility

O parâmetro `walltime` agora é **scheduler-agnostic**:

| Scheduler | Formato  | Exemplo       |
|-----------|----------|---------------|
| SLURM     | HH:MM:SS | `2:00:00`     |
| PBS/Torque| HH:MM:SS | `02:00:00`    |
| Genérico  | HH:MM:SS | `1:30:45`     |

O xespresso passa o valor diretamente para o scheduler, que o interpreta apropriadamente.

## Próximos Passos

1. **Executar teste periodicamente**: Garante que mudanças futuras não quebrem jobfile generation
2. **Adicionar a CI/CD**: Validar em cada commit
3. **Testar com jobs reais**: Verificar que jobs são submetidos e executados corretamente
4. **Documentar formatos de walltime**: Adicionar guia para cada scheduler específico
