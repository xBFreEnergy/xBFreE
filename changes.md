
# TODO
Cambios estructurales
- [ ] el rewrite-output cambió por la estructura de carpetas
  - [ ] es necesario arreglar esto
  - [ ] hacerlo compatible con las últimas versiones de gmx_MMPBSA
  - [ ] 



Cambios añadidos
- [ ] añadir una nueva capa al archivo de resultados binario
  - [ ] añadir el id del sistema usando hash a partir de las opciones -cs y -cp?
  - [ ] añadir la versión del programa. Detectar si es una rewrite-output de gmx_MMPBSA
  - [ ] añadir el método de cálculo (mmpbsa, lie...)
  - [ ] comprobar si el binario termina en mmxsa o xbfree
- [ ] 


API
- [ ] load_file
  - [ ] debe aceptar archivos individuales o listas
  - [ ] debe asignar el tipo de cálculo al sistema si está anotado
- [ ] el API debe almacenar los CalculationsSystem?

