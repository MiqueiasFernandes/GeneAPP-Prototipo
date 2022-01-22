<template>
  <Tabs>
    <template #tab1-title>Junction reads coverage</template>
    <template #tab1>
      <ul>
        <li>
          grafido de barras por amostra (juntar replicas) do tamanho da juncao
          amostrado empilhado por tipo basico de evento
        </li>
        <li>
          grafico de linha do tp basico e evt com a profundidade de reads
          down/up stream do sitio alternativo por amostra
        </li>
      </ul>
    </template>
    <template #tab2-title> Events overview </template>
    <template #tab2>
      <ul>
        <li>barras por tipo de evento</li>
        <li>volcano de evento, signifig por ∂psi</li>
      </ul>
    </template>

    <template #tab3-title> IsoStar </template>
    <template #tab3>
      <ul>
        <li>grafico de genes por amostra por iso com mais tpm</li>
      </ul>
    </template>

    <template #tab4-title> Tabela </template>
    <template #tab4>
      <ul>
        <li>tabela dos eventos chr,gene,iso,grupo,tipo,evento,psi,fdr,posa,posb,iso2,tpmgene1,tpmgene2,tpmiso1,tpmiso2,tpmiso21, tpmiso22,lengene, lenexon,lencds,lenintro,qtdmrna,qtdex,qtdin,emcds,cdsaltsize</li>
        <li>filtrar colunas</li>
        <li>filtrar por carac</li>
      </ul>
    </template>

    <template #tab5-title> samples </template>
    <template #tab5>
      <ul>
        <li>violin do contexto genomico das features total, por amostra, por tipo evento</li>
      </ul>
    </template>
  </Tabs>
</template>

<script>
export default {};
</script>

<style>
</style>