# TFM-CarlosYanesPerez
Códigos implementados para el TFM de Carlos Yanes Pérez, con título: Ataques cuánticos contra criptografía resistente a la cuántica.

En el directorio se encuentran varios notebooks diseñados para ser ejecutados en Google Colab.

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/)

1. S-AES
   * Ataque VQA. Implementación del ataque VQE al S-AES basado en el texto de Wang, Z., Zheng, M., Wu, J. et al. *Reducing quantum resources for attacking S-AES on quantum devices.* npj Quantum Inf 11, 157 (2025). (DOI: [10.1038/s41534-025-01106-w](https://www.nature.com/articles/s41534-025-01106-w)). Además, se ha hecho uso de las notas suplementarias de dicha publicación.
   * Ataque mediante Grover. Se implementa el S-AES cuántico descrito en Escanez-Exposito & Caballero-Gil (2026). *Quantum simulation of Boolean logic, digital circuits, and cryptographic schemes*. Journal of Logic and Computation. DOI: [10.1093/jigpal/jzaf019](https://academic.oup.com/jigpal/article/34/1/jzaf019/8443240). Posteriormente, de implementa el algoritmo de Gover contra parte de la clave buscada.

2. LWE
   * Ataque HAWI. Implementación del algoritmo HAWI propuesto en Zheng et al. (2025). *Quantum-classical hybrid algorithm for solving the learning-with-errors problem on NISQ devices*. Communications Physics. DOI: [10.1038/s42005-025-02126-w](https://doi.org/10.1038/s42005-025-02126-w).
   * Ataque mediante *Quantum Sieving*. Implementación del algoritmo de cribado cuántico propuesto en: Rojas et al. (2025). *Quantum Sieving to solve the LWE Problem*. ITEFI-CSIC, Madrid.

3. Cifrado Afín
   * Ataque mediante Grover. Inspirado por el trabajo de M. M. Mathews, P. V and V. Ajith, "Quantum Cryptanalysis of Affine Cipher," in IEEE Journal on Emerging and Selected Topics in Circuits and Systems, vol. 14, no. 3, pp. 507-519, Sept. 2024, doi: [10.1109/JETCAS.2024.3428436](https://doi.org/10.1109/JETCAS.2024.3428436).

4. Implementacion del SAES

5. Modos de cifrado AES
