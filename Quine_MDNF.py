# 10.04.2017 by galaxid3d
# Находит все МКНФ/МДНФ заданной ф-ии методом Квайна. Ф-ию можно вводить как строку из ('1' и '0' и '*') или из (чисел; если надо, то число*)

def count_elements(i):  # определяет необходимое кол-во элементов для выражения ф-ии
    if i < 0: i = -i
    if i < 2: return 1
    result = 0
    k = 1
    while k <= i:
        result += 1
        k = k << 1
    return result

def isValid_length(length):  # ввели ли допустимую длину ф-ции
    one_count = 0  # количество едениц в числе (должна быть одна - тогда это степень двойки)
    s = str(bin(length))[2:]  # убираем "0b" из числа
    for i in s:
        if i == "1":
            one_count += 1
            if one_count > 1: return False
    return one_count

def make_term(k, N, isDNF=True):  # делает терм из целого числа заданной длинны
    result = []
    s = str(bin(k))[2:]  # убираем "0b" из числа
    for i in range(N - len(s)):
        if isDNF: result.append(-(1 + i))
        else: result.append(1 + i)
    for i in range(len(s), 0, -1):
        if s[len(s) - i] == str(int(isDNF)):
            result.append(N + 1 - i)
        else:
            result.append(-(N + 1 - i))
    return result

def make_join(term1, term2):  # делает склейку: xy + xyz = xy
    if len(term1) != len(term2): return False
    isDimension = False
    result = []
    for i in range(len(term1)):
        if abs(term1[i]) == abs(term2[i]):
            if term1[i] == term2[i]:
                result.append(term1[i])
            else:
                if isDimension:
                    return False
                else:
                    isDimension = True
        else:
            return False  # если разные числа
    if isDimension:
        return result
    else:
        return False

def make_merge(term1, term2):  # делает поглощение: xy + x!y = x
    if abs(len(term1) - len(term2)) != 1: return False
    if len(term1) < len(term2):
        tmp = term1
        term1 = term2
        term2 = tmp
    k = 0
    for j in range(len(term1)):
        for i in range(len(term2)):
            if term1[j] == term2[i]: k += 1
    if k == len(term2):
        return term2
    else:
        return False

def print_MDNF(MDNF, isDNF=True):  # печатает ДНФ/КНФ
    for j in range(len(MDNF)):
        if not isDNF: print(end="(")
        for i in range(len(MDNF[j])):
            if MDNF[j][i] < 0: print(end="!")
            print(end="x" + str(abs(MDNF[j][i])))
            if (not isDNF) and (i + 1 < len(MDNF[j])): print(end=" √ ")
        if not isDNF: print(end=")")
        if j + 1 < len(MDNF):
            if isDNF:
                print(end=" √ ")
            else:
                print(end=" & ")
    print()

def isIn_term(term1, term2):  # входит ли один терм в другой
    if len(term1) < len(term2):
        tmp = term1
        term1 = term2
        term2 = tmp
    for element in term2:
        if term1.count(element) == 0: return False
    return True

def make_all_joints(super_DNF):  # делает все склейки
    result = False
    out_DNF = []
    used = [False for _ in range(len(super_DNF))]
    for j in range(len(super_DNF)):
        for i in range(j + 1, len(super_DNF)):
            tmp_Term = make_join(super_DNF[j], super_DNF[i])
            if tmp_Term:
                if not (tmp_Term in out_DNF):
                    out_DNF.append(tmp_Term)
                result = used[j] = used[i] = True
    if result:
        for i in range(len(used)):
            if not used[i]:
                out_DNF.append(super_DNF[i])
        return out_DNF
    else:
        return False

def make_all_merges(super_DNF):  # делает все поглощения
    result = False
    out_DNF = []
    used = [False for _ in range(len(super_DNF))]
    for j in range(len(super_DNF)):
        for i in range(j + 1, len(super_DNF)):
            tmp_Term = make_merge(super_DNF[j], super_DNF[i])
            if tmp_Term:
                if not (tmp_Term in out_DNF):
                    out_DNF.append(tmp_Term)
                result = used[j] = used[i] = True
    if result:
        for i in range(len(used)):
            if not used[i]:
                out_DNF.append(super_DNF[i])
        return out_DNF
    else:
        return False

def copy_MDNF(from_MDNF):  # копирует из одной в другую
    return [term for term in from_MDNF]

def get_var_count(DNF):  # считает кол-во букв в ДНФ
    result = 0
    for Term in DNF: result += len(Term)
    return result

def record_update(MDNF):  # обновляем лучшую МДНФ
    MDNF_len = get_var_count(MDNF)
    global record_len
    global record
    if MDNF_len == record_len:
        record.append(copy_MDNF(MDNF))
    if MDNF_len < record_len:
        record_len = MDNF_len
        record = [copy_MDNF(MDNF)]

def isAll_used(MDNF):  # покрывает ли МДНФ м-цу импликативности
    for i in range(len(super_DNF)):
        inColumn = False
        for j in range(len(MDNF)):
            inColumn = inColumn or isIn_term(super_DNF[i], MDNF[j])
            if inColumn: break
        if not inColumn: return False
    return inColumn

def find_MDNF(MDNF, pos, posT, count):  # вычисление МДНФ: пробует все термы в разном кол-ве. Поэтому матрица импликативности не нужна
    if pos == count:
        if isAll_used(MDNF): record_update(MDNF)
        return
    for i in range(posT, len(short_DNF)):
        MDNF[pos] = short_DNF[i]
        find_MDNF(MDNF, pos + 1, i + 1, count)

def get_all_DNF(poses, pos, posT, count, isDNF):  # если ввели ДНФ со звёздочками, то она находит всевозможные ДНФ из них
    global orig_func
    if pos == count:
        tmp_s = orig_func
        for posN in poses:
            tmp_s = tmp_s[:posN] + isDNF + tmp_s[posN + 1:]
        all_DNF.append(tmp_s)
        return
    for i in range(posT, len(pos_stars)):
        poses[pos] = pos_stars[i]
        get_all_DNF(poses, pos + 1, i + 1, count, isDNF)

# Input function
s = input("Is DNF: ")
isDNF = s != "0" and s != "False" and s != ""
s = input("Input as numbric: ")
if s != "0" and s != "False" and s != "":
    print(end="Input function: ")
    sFunc = list(map(str, input().split()))  # или вводим только нужные числа
    func = [int(c.replace("*", "")) for c in sFunc]
    s = ""
    for j in range(2 ** count_elements(max(int(input("Input max element: ")), max(func)))):  # если ошиблись и ввели меньше, чем есть
        for i in range(len(func)):
            if func[i] == j:
                if "*" in sFunc[i]: s += "*"
                else: s += '1'
                break
        else: s += '0'
else:
    s = input("Input function: ")  # вводим саму ф-ю: 0110

# Start Quine algorithm
if isValid_length(len(s)):
    # Get All DNF by stars: 0**1=[0011,0101,0111/0001-в зависимости от isDNF]
    pos_stars = []
    for i in range(len(s)):
        if s[i] == "*": pos_stars.append(i)
    orig_func = s.replace("*", str(int(not isDNF)))
    all_DNF = [orig_func]  # размер(т.е. всевозможные комбинации из '*') = 2**(кол-во '*')
    for count in range(1, len(pos_stars) + 1):
        poses = [0 for i in range(count)]
        get_all_DNF(poses, 0, 0, count, str(int(isDNF)))
    N = count_elements(len(s) - 1)
    best_MDNF = []  # если ввели ф-ию с "*", тут хранятся все ф-ии с наим. кол-вом и эл-тов и их МДНФ
    best_len = len(s) ** 2
    funcs = []  # лучшие ф-ии (с минимальной МДНФ)
    for s in all_DNF:
        # Make super_DNF - совершенная ДНФ
        super_DNF = []
        for i in range(len(s)):
            if s[i] == str(int(isDNF)):
                super_DNF.append(make_term(i, N, isDNF))

        # Make short_DNF (Joins and Mergings) - сокращённая ДНФ
        short_DNF = [Term for Term in super_DNF]
        isChange = True
        while isChange:
            isChange = False
            tmp_DNF = make_all_joints(short_DNF)  # временная ДНФ, куда запишется рез-тат склейки/поглощения
            if tmp_DNF:
                isChange = True
                short_DNF = copy_MDNF(tmp_DNF)
            tmp_DNF = make_all_merges(short_DNF)
            if tmp_DNF:
                isChange = True
                short_DNF = copy_MDNF(tmp_DNF)

        # Make MDNF
        record = [[Term for Term in short_DNF]]  # создаём список самых коротких (лучших) ДНФ
        record_len = get_var_count(record[0])
        for count in range(1, len(short_DNF)):
            MDNF = [[] for i in range(count)]
            find_MDNF(MDNF, 0, 0, count)

        # Update best_MDNF
        MDNF_len = get_var_count(record[0])
        if MDNF_len == best_len:
            best_MDNF.append([copy_MDNF(MDNF) for MDNF in record])
            funcs.append(s)
        if MDNF_len < best_len:
            best_len = MDNF_len
            best_MDNF = [[copy_MDNF(MDNF) for MDNF in record]]
            funcs = [s]

    # Print All MDNF
    for j in range(len(best_MDNF)):
        s = "Function: " + funcs[j]
        print("_" * len(s));
        print()
        print(s)
        for MDNF in best_MDNF[j]:
            if isDNF:
                print(end="MDNF: ")
            else:
                print(end="MKNF: ")
            print_MDNF(MDNF, isDNF)
else:
    print("Incorrect length of function!")