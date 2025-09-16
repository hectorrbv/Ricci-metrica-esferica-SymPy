# Ricci en métrica esféricamente simétrica (SymPy)

Calculadora simbólica del **tensor de Ricci** $R_{\mu\nu}$ para la métrica esféricamente simétrica

$$
ds^2=-B(r)\,dt^2 + A(r)\,dr^2 + r^2\,d\theta^2 + r^2\sin^2\theta\,d\phi^2,
$$

implementada en **SymPy**. El script imprime:

* Los **símbolos de Christoffel** usados (uno por uno).
* Las **componentes de $R_{\mu\nu}$** calculadas desde la **definición** (sumas con índices).
* Un bloque **LaTeX** con las fórmulas cerradas esperadas para:

  $$
  \begin{aligned}
  R_{tt} &= -\frac{B''}{2A} - \frac{B'}{4A}\!\left(\frac{B'}{B} - \frac{A'}{A}\right) + \frac{1}{r}\frac{B'}{A},\\
  R_{rr} &= \phantom{-}\frac{B''}{2B} + \frac{B'}{4B}\!\left(\frac{B'}{B} - \frac{A'}{A}\right) + \frac{1}{r}\frac{A'}{A},\\
  R_{\theta\theta} &= 1 - \frac{1}{A} + \frac{r}{2A}\!\left(\frac{A'}{A} - \frac{B'}{B}\right),\qquad
  R_{\phi\phi}=\sin^2\theta\,R_{\theta\theta}.
  \end{aligned}
  $$

> ⚠️ **Convenciones**: El resultado depende de la **definición de Riemann/Ricci**. Este repo usa la convención
> $R_{\mu\nu}=\partial_\nu \Gamma^\rho_{\mu\rho}-\partial_\rho \Gamma^\rho_{\mu\nu}+\Gamma^\rho_{\nu\sigma}\Gamma^\sigma_{\mu\rho}-\Gamma^\rho_{\rho\sigma}\Gamma^\sigma_{\mu\nu}$
> (típica en varios textos). Cambiar la convención invierte signos en algunos términos.

---

## Requisitos

* Python 3.9+
* [SymPy](https://www.sympy.org/) (>= 1.12)
* NumPy (opcional; el script lo importa pero no es estrictamente necesario)

Instalación rápida:

```bash
pip install sympy numpy
```

---

## Uso

Ejecuta el script directamente:

```bash
python ricci_esferico.py
```

Salida esperada (extracto):

* Listado de $\Gamma^\mu{}_{\nu\rho}$ no nulos.
* Cálculo de $R_{\mu\nu}$ “desde la definición”.
* **Bloque LaTeX** con las componentes cerradas:

  ```latex
  \begin{align}
  R_{tt} &= \cdots \\
  R_{rr} &= \cdots \\
  R_{\theta\theta} &= \cdots \\
  R_{\phi\phi} &= \cdots
  \end{align}
  ```
* Matriz diagonal de $R_{\mu\nu}$ (en el sistema de coordenadas dado).

> 👉 Copia/pega el bloque LaTeX directamente a tu nota/reporte.

---

## Cómo cambiar la métrica (A y B)

Por defecto $A=A(r)$, $B=B(r)$ son funciones simbólicas. Para probar casos concretos, puedes **sustituir** al final del script:

```python
M = sp.symbols('M', positive=True)
subs_dict = {
    A: 1/(1-2*M/r),   # Ej. Schwarzschild: A=1/f
    B: (1-2*M/r)      #                     B=f
}
print(sp.simplify(R11.subs(subs_dict)))
```

> **Unidades**: usualmente $G=c=1$.

---

## Verificación rápida

* La salida “directa” usa las fórmulas cerradas.
* La salida “desde definición” usa la suma explícita sobre $\rho,\sigma$ y derivadas parciales en $r$ y $\theta$.
* Deben coincidir simbólicamente para $R_{tt},R_{rr},R_{\theta\theta}$ y $R_{\phi\phi}$ (este último igual a $\sin^2\theta\,R_{\theta\theta}$).

Si cambias la **convención** de curvatura (o el signo en la métrica), es normal observar **signos distintos** en términos como $B''$ o los cruces $A'B'$.

---

## Extensiones útiles (opcional)

Puedes añadir al final del script:

**Escalar de Ricci $R$:**

```python
gtt, grr = -B, A
gthth, gphph = r**2, r**2*sp.sin(theta)**2
g_inv = sp.diag(-1/B, 1/A, 1/r**2, 1/(r**2*sp.sin(theta)**2))
R_scalar = sp.simplify(g_inv[0,0]*R00 + g_inv[1,1]*R11 + g_inv[2,2]*R22 + g_inv[3,3]*R33)
print("R =", sp.latex(R_scalar))
```

**Tensor de Einstein $G_{\mu\nu} = R_{\mu\nu} - \tfrac12 g_{\mu\nu} R$:**

```python
g = sp.diag(-B, A, r**2, r**2*sp.sin(theta)**2)
G = sp.simplify(sp.Matrix([
    [R00 - sp.Rational(1,2)*g[0,0]*R_scalar, 0, 0, 0],
    [0,   R11 - sp.Rational(1,2)*g[1,1]*R_scalar, 0, 0],
    [0, 0, R22 - sp.Rational(1,2)*g[2,2]*R_scalar, 0],
    [0, 0, 0, R33 - sp.Rational(1,2)*g[3,3]*R_scalar],
]))
print("G_{\\mu\\nu} =", sp.latex(G))
```

---

## Problemas frecuentes

* **“No coincide con mis apuntes”** → Revisa la **convención** de $R^\rho{}_{\sigma\mu\nu}$ usada en tu curso vs. la de este script. Cambiarla invierte signos en varios términos.
* **“No aparecen términos $1/r$”** → Asegúrate de no eliminar los $\Gamma$ angulares $\Gamma^\theta_{r\theta}=\Gamma^\phi_{r\phi}=1/r$, $\Gamma^r_{\theta\theta}=-r/A$, $\Gamma^r_{\phi\phi}=-r\sin^2\theta/A$.
* **“SymPy tarda/satura”** → Usa `simplify` solo en la salida final, evita simplificar cada término intermedio.

---

## Licencia

MIT. Úsalo libremente con atribución.

---

## Cita

Si te sirve en un trabajo o repo, puedes citar así:

> *“Ricci en métrica esféricamente simétrica (SymPy)”, Héctor Becerril (2025), GitHub repository.*

