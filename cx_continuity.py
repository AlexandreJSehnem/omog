# Aluno: Alexandre Jeronimo Sehnem

### COMO EXECUTAR O CODIGO ###
#   Requer:
#       - Python3 instalado;
#       - pip instalado (gerenciador de pacotes do python);
#       - biblioteca numpy, comando: pip install numpy
#       - biblioteca matplotlib, comando: pip install matplotlib
#
#   Comando para executar: python cx_continuity.py


# Descrição: este código implementa uma curva B-sppline de grau 6, uma bezier de grau 3 e então unifica as duas com continuidade C0


import numpy as np
import matplotlib.pyplot as plt

pontos = np.array([[1, 1], [2, 4], [4, 5], [6, 3], [7, 0], [9, 2], [10, 5], [11, 4]])

vetor_nos = np.array([0, 0, 0, 0, 0, 0, 0, 1, 2, 2, 2, 2, 2, 2, 2])
grau = 6

# algoritmo de De Boor, faz o calculo da funcao para encontrar cada um dos pontos da Spline
def de_boor(i, grau, vetor_nos, t):
    if grau == 0:
        return 1.0 if vetor_nos[i] <= t < vetor_nos[i + 1] else 0.0

    denom1 = vetor_nos[i + grau] - vetor_nos[i]
    denom2 = vetor_nos[i + grau + 1] - vetor_nos[i + 1]

    sum1 = 0.0 if denom1 == 0.0 else (t - vetor_nos[i]) / denom1 * de_boor(i, grau - 1, vetor_nos, t)
    sum2 = 0.0 if denom2 == 0.0 else (vetor_nos[i + grau + 1] - t) / denom2 * de_boor(i + 1, grau - 1, vetor_nos, t)

    return sum1 + sum2

# calcula os pontos da curva
def bspline(pontos, vetor_nos, grau, t):
    n = len(pontos) - 1
    pontos_curva = np.zeros_like(pontos[0])

    for i in range(n + 1):
        funcao_b = de_boor(i, grau, vetor_nos, t)
        pontos_curva = pontos_curva + (funcao_b * pontos[i])

    return pontos_curva

# cria um vetor de numeros tal que, para todo n pertencente a t_values e vetor_nos[grau] <= n <= vetor_nos[-(grau+1)]
# alem disso, cada valor de t_values esta a uma distancia fixa um do outro
vetor_t = np.linspace(vetor_nos[grau], vetor_nos[-(grau + 1)], num=100)

# utimo ponto eh aproximado
vetor_t[-1] = vetor_t[-1] - 0.00001

curva_final = [bspline(pontos, vetor_nos, grau, t) for t in vetor_t]

# pontos b-spline
x_spline = [ponto[0] for ponto in curva_final]
y_spline = [ponto[1] for ponto in curva_final]

############## COMEÇO BEZIER ###############

# Algoritmo que aplica a função para calcular os pontos da bezier de grau 3
def bezier(P0, P1, P2, P3, num_points):
    t = np.linspace(0, 1, num_points)
    x = (1 - t)**3 * P0[0] + 3 * (1 - t)**2 * t * P1[0] + 3 * (1 - t) * t**2 * P2[0] + t**3 * P3[0]
    y = (1 - t)**3 * P0[1] + 3 * (1 - t)**2 * t * P1[1] + 3 * (1 - t) * t**2 * P2[1] + t**3 * P3[1]
    return x, y

# Pontos de Controle
P0 = [1, 3]
P1 = [3, 3]
P2 = [4, 4]
P3 = [8, 4]

# Qtd. Pontos da Curva
num_points = 100

# Gera Bezier
x_bezier, y_bezier = bezier(P0, P1, P2, P3, num_points)

# Função de plot das duas curvas juntas
def plot_curvas(title):
	plt.plot(x_spline, y_spline, label='B-spline')
	plt.plot(pontos[:, 0], pontos[:, 1], marker='o', color='r', label='Pontos de Controle')

	plt.plot(x_bezier, y_bezier, label="Bezier")
	plt.plot([P0[0], P1[0], P2[0], P3[0]], [P0[1], P1[1], P2[1], P3[1]], marker='o', c="g", label="Pontos de Controle")

	plt.legend()
	plt.xlabel("X")
	plt.ylabel("Y")
	plt.title(title)
	plt.grid(True)
	plt.show()

# plot antes de implementar C0
plot_curvas("Plot Inicial")



############ IMPLEMENTANDO C0 ############
# calculando delta para transladar os pontos de controle da curva
x_delta = pontos[-1][0] - P0[0]
y_delta = pontos[-1][1] - P0[1]

# transladando os pontos
P0[0] = P0[0] + x_delta
P0[1] = P0[1] + y_delta
P1[0] = P1[0] + x_delta
P1[1] = P1[1] + y_delta
P2[0] = P2[0] + x_delta
P2[1] = P2[1] + y_delta
P3[0] = P3[0] + x_delta
P3[1] = P3[1] + y_delta

# recriando a curva de bezier com os novos pontos de controle
x_bezier, y_bezier = bezier(P0, P1, P2, P3, num_points)

# plot de c0
plot_curvas("Plot C0")


# inicio continuidade c1 e c2
#calcula derivada d
def derivada_spline(i, grau, vetor_nos, t, d):
    if grau == 0:
        return 1.0 if vetor_nos[i] <= t < vetor_nos[i + 1] else 0.0
    
    denom1 = vetor_nos[i + grau - 1] - vetor_nos[i]
    denom2 = vetor_nos[i + grau] - vetor_nos[i+1]

    if d > 0:
        sum1 = 0.0 if denom1 == 0.0 else (grau - 1) / denom1 * derivada_spline(i, grau - 1, vetor_nos, t, d-1)
        sum2 = 0.0 if denom2 == 0.0 else (grau - 1) / denom2 * derivada_spline(i + 1, grau - 1, vetor_nos, t, d-1)
    else:
        return de_boor(i, grau, vetor_nos, t)

    return sum1 - sum2


# calcula os pontos da curva com derivada d
def bspline_pontos_derivada(pontos, vetor_nos, grau, t, d):
    n = len(pontos) - 1
    pontos_curva = np.zeros_like(pontos[0])

    for i in range(n + 1):
        funcao_b = derivada_spline(i, grau, vetor_nos, t, d)
        pontos_curva = pontos_curva + (funcao_b * pontos[i])

    return pontos_curva

# faz união com continuidade c1:
def unify_c1():
    global pontos, vetor_nos, grau, vetor_t

    derivada_spline = bspline_pontos_derivada(pontos, vetor_nos, grau, (vetor_t[-1]-0.00001), 1)

    P1[0] = (derivada_spline[0]/(4-1)) + P0[0]
    P1[1] = (derivada_spline[1]/(4-1)) + P0[1]

    print(derivada_spline)

unify_c1()
x_bezier, y_bezier = bezier(P0, P1, P2, P3, num_points)
plot_curvas("c1")


# faz união com continuidade c2:
def unify_c2():
    global pontos, vetor_nos, grau, vetor_t

    derivada_spline = bspline_pontos_derivada(pontos, vetor_nos, grau, (vetor_t[-1]-0.00001), 2)

    P2[0] = (derivada_spline[0]/(4*(4-1))) + (2*P1[0] - P0[0])
    P2[1] = (derivada_spline[1]/(4*(4-1))) + (2*P1[1] - P0[1])

    print(derivada_spline)

unify_c2()
x_bezier, y_bezier = bezier(P0, P1, P2, P3, num_points)
plot_curvas("c2")


