import subprocess

def compilar_e_executar_cpp():
    # Compilação
    subprocess.run(['g++', 'Algoritimo_Genetico_1.cpp', '-o', 'meu_programa'])

    # Execução
    subprocess.run(['./meu_programa'])

compilar_e_executar_cpp()
