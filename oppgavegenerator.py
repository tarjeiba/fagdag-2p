import math
import pdb
 
import random as rnd
import numpy as np
import matplotlib.pyplot as plt
import statistics as stat

from matplotlib.ticker import AutoMinorLocator

# TODOS:
# - TODO :: Diagramtolkning: beholder a-versjonen, men gjør slik at den ene grafen får litt større spredning. 
# - TODO :: Endre fasitgenerering til å ta høyde for at frekvenstabell- og analyseoppgavene er endret
# - TODO :: plot_list har en ekstremt hackete løsning for å endre spredningen. Denne bør egentlig ikke bruke
#           histogram for å finne data

colors = {
    'blue' : (0.4, 0.5977, 1.0),
    'pink' : (0.8, 0.7, 1.0),
    }
 
groups = ['test']
levels = ['a']
number_of_students = 2
 
replacements = {}
 
words_easy = ['median', 'typetall', 'gjennomsnitt', 'kvartilbredde']
words_hard = ['statistisk modell', 'historgram', 'grupperte data', 'varians', 'standardavvik']
 
 
template_file_name = 'oppgave-mal.org'
 
template_file = open(template_file_name, 'r', encoding="utf-8")
template = template_file.read()
template_file.close()
 
 
def median_list(lst):
    """Return the median value of a list."""
    srtlst = sorted(lst)
 
    length = len(srtlst)
 
    if length % 2:
        return srtlst[math.floor(len(srtlst) / 2)]
    else:
        middle = length/2
         
        return sum(srtlst[(length//2-1):(length//2+1)])/2

def range_lst(lst):
    """Return the range of a list of numbers."""
    return max(lst) - min(lst)

def interquartile_range_list(lst):
    """Return the interquartile range of a list of numbers."""
    lower = np.percentile(lst, 25)
    higher = np.percentile(lst, 75)
    return higher - lower
 
def range_freq_table(lst_of_tups):
    """Return the range of the values in a frequency table."""
    values = [tup[0] for tup in lst_of_tups]
    counts = [tup[1] for tup in lst_of_tups]
    low_ind = 0
    high_ind = -1
    while True:
        if counts[low_ind] != 0:
            lowest = values[low_ind]
            break
        low_ind += 1
    while True:
        if counts[high_ind] != 0:
            highest = values[high_ind]
            break
        high_ind -= 1
 
    return highest - lowest
 
 
def median_freq_table(lst_of_tups):
    """Return the median of the values in a frequency table."""
    lst = []
    for tup in lst_of_tups:
        lst.extend([tup[0]] * tup[1])
    return median_list(lst)
 
 
def mode_freq_table(lst_of_tups):
    """Return the mode of the values in in a frequency table."""
    values = [tup[0] for tup in lst_of_tups]
    counts = [tup[1] for tup in lst_of_tups]
    max_count = max(counts)
    return values[counts.index(max_count)]
 
def create_line_data(n, stddev=1000, sol=False, dictkey=None):
    """Return a list of n (year, value) tuples."""
 
    start_year = 1996
    years = range(start_year, start_year + n)
 
    coefficient = rnd.randint(1, 100)
 
    if rnd.random() > 0.5:
        coefficient *= -1
 
    k_value = rnd.random()
    if rnd.random() > 0.5:
        k_value *= -1
 
    values = [k_value * t + coefficient for t in range(n)]
    values = [i + rnd.gauss(0, stddev) for i in values]
 
    return list(zip(years, values))
 
def average_list(lst):
    """Return average of all elements in list lst."""
    return sum(lst)/len(lst)
 
def average_freq_dice(lst):
    """Return average of all elements in list of dice throws
    where index is value - 1.
    """
    sums = 0
    tot = 0
    for val, freq in enumerate(lst):
        sums += (val + 1) * freq
        tot += freq
 
    return sums/tot
 
def random_list_easy(n, maximum=6, sol=False, dictkey=None):
    """Return a list of length n where the numbers are "easy" to use for
    statistical analysis."""
    value = [rnd.randint(1, maximum) for _ in range(n)]
 
    if sol:
        avg = average_freq_dice(value)
        repl[dictkey] = ("Gjennomsnittet er {0}."
                                 "".format(avg))
    return list(map(str, value))
 
def random_list_hard(n=32, mu=15, sigma=10, sol=False, dictkey=None):
    """Return a list of length n where the value of each member is
    Gaussian distributed with an expected value of mu and a std.dev. of sigma.
 
    n has a default value of 32 due to it being a typical class size.
    """
    values = [rnd.gauss(mu, sigma) for _ in range (n)]
    if sol:
        repl[dictkey] = ("Medianen er {:.2f}.\n"
                         "Gjennomsnittet er {:.2f}.\n"
                         "Kvartilbredden er {:.2f}.\n"
                         "Standardavviket er {:.2f}.\n"
                         .format(
                             median_list(values),
                             average_list(values),
                             interquartile_range_list(values),
                             np.std(values, ddof=1) # TODO :: ddof=? Spør
                         ))
    return ["{0:.2f}".format(val) for val in values]
 
 
def histogram_addition():
    return "- Ut fra verdiene i histogrammet over, hva vil du tro mengden pastiller i krukka var?\n"
 
 
def histogram_data(num=25, level="test", seed_string="test", sol=False, dictkey=None):
    # if sol:
    #     repl[dictkey] = "Fasiten er {}.".format(num)
    return ([rnd.gauss(25, 10) for _ in range(num)], [0, 5, 25, 30, 40, 50])
 
 
# TODO: Måten histogram_data blir kalt på, gjør at den i praksis kun setter 
 
def plot_histogram(data, file_name="test_hist.png", seed_string="test", sol=False, dictkey=None):
    """Create a histogram.
     
    data is a tuple of two lists:
    g_list = [a, b, c, d, e, ...] # length n + 1
    f_list = [f_ab, f_bc, f_cd, ...] # length n
    """
 
    num = rnd.randint(25, 50)
    figure_name = file_name[:-4] + "_hist.png"
    plot_figure = plt.figure()
 
    ax = plot_figure.add_subplot(111)
    hist, bins_ = np.histogram(data[0], bins=data[1], density=True)
    hist_data = [num * x for x in hist]
    widths = [bins_[i+1] - bins_[i] for i in range(len(bins_)-1)]
    lefts = [bins_[i] + widths[i]/2 for i in range(len(widths))]
 
    ax.bar(lefts, hist_data, width=widths)
 
    for x, y, w in zip(lefts, hist_data, widths):
        ax.text(x, y, "{:.3f}".format(y))
 
    ax.grid(linestyle="dashed", which='both')
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.set_xticks(bins_)
    ax.set_xlabel("Antall pastiller")
    ax.set_ylabel("Høyde")
    plot_figure.savefig("./figurer/" + figure_name)
    plt.close(plot_figure)
 
    if sol:
        repl[dictkey] = "Fasiten er {}.".format(num)
 
    return insert_image(figure_name)
 
 
def plot_list(name="test.png", start=None, end=None, stepsize=None):
    # data1, bins = np.histogram([rnd.triangular(0.0, 10, mode=4) for _ in range(1000)], 10)
    # data2, _ = np.histogram([rnd.triangular(0.0, 10, mode=4) for _ in range(4000)], 10)
    err = 8
    data1, bins = np.histogram([rnd.gauss(0.0, 1) for _ in range(1000)], 10)
    data2, _ = np.histogram([rnd.gauss(0.0, 1) for _ in range(4000)], 10 + 2 * err)

    plot_list_bar(data1, data2[err // 2  + 2: 12 + err // 2], name=name, start=0, end=10, stepsize=stepsize)
 
def plot_list_line(lst1, lst2=None, name='test_line.png'):
    """Create a line graph of lst1 (and or 2)."""
 
    plot_figure = plt.figure()
    ax = plot_figure.add_subplot(111)
    years = [i[0] for i in lst1]
    data = [i[1] for i in lst1]
    ax.plot(years, data, color=colors['blue'])
    if lst2:
        years2 = [i[0] for i in lst2]
        data2 = [i[1] for i in lst2]
        ax.plot(years2, data2, color=colors['pink'])
 
    ax.set_xticks(years[::5])
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    plot_figure.savefig("./figurer/" + name)
    plt.close(plot_figure)
    return None
 
 
def plot_list_bar(lst1, lst2=None, name='test.png', start=None, end=None, stepsize=None, xticks=None):
    """Create a bar graph of lst (and lst2). The function returns void,
    creates a file.
    """
     
    plot_figure = plt.figure()
    ax = plot_figure.add_subplot(111)
    n = len(lst1)
    ind = np.arange(n)
    width = 0.35
    d = 0.03
 
    if lst2 is not None:
        ax.bar(ind - width/2 - d, lst1, width, color=colors['blue'])
        ax.bar(ind + width/2 + d, lst2, width, color=colors['pink'])
    else:
        ax.bar(ind, lst1, width, color=colors['blue'])
 
    if stepsize and not xticks:
        ax.xaxis.set_ticks(np.arange(start, end, stepsize))
 
    if xticks is not None:
        ax.xaxis.set_ticks(xticks)
         
    plot_figure.savefig("./figurer/" + name)
    plt.close(plot_figure)
    return None
 
def insert_image(filename):
    """Inserts an org-link to an image file."""
    return  "[[./figurer/" + filename + "]]"
 
 
def draw_words(word_list, n=3):
    """Replace value in replacements_dict at replacement_key with n random values from word_list."""
    drawn_words = rnd.sample(word_list, n)
    drawn_string = ", ".join(drawn_words[:-1])
    drawn_string += " og " + drawn_words[-1]
    return drawn_string
 
def populate_freq_table_easy(n, t=6, seed_string="test", dictkey=None):
    """Populate a frequency table with n dice throws from 1 to 6."""
    throws = [0] * t
    return_string = ""
    rnd.seed(seed_string)
 
    for i in range(n):
        throws[rnd.randint(0, t-1)] += 1
    for l in range(t):
        return_string += "| {0} | {1} |\n".format(l+1, throws[l])
 
    freq_tup = list(zip(range(1, t+1), throws))
    avg = average_freq_dice(throws)
    med = median_freq_table(freq_tup)
    typ = mode_freq_table(freq_tup)
    var = range_freq_table(freq_tup)
    solutions[dictkey] = ("Gjennomsnittet er {0:.2f}.\n"
                         "Medianen er {1}.\n"
                         "Typetallet er {2}.\n"
                         "Variasjonsbredden er {3}.\n"
                         "".format(avg, med, typ, var))
    return return_string
 
 
def populate_freq_table_hard(n=None, t=8, seed_string="test", dictkey=None):
    """Populate a frequency table with scores from 0 to 100.
    """
    bin_levels = [0, 20, 40, 60, 80, 95, 100]
    bins = ['[{0}, {1}>'.format(bin_levels[i], bin_levels[i+1]) for i in range(len(bin_levels)-1)]
    return_string = ''
    hist = []
 
    for binn in bins:
        ran = rnd.randint(0, 10)
        hist.append(ran)
        return_string += "| {0} | {1} |\n".format(binn, ran)
     

    # Creating solutions
        
    center_list = centers(bin_levels)
    linspace_list = []
    for i in range(len(bin_levels) - 2):
        temp = np.linspace(bin_levels[i], bin_levels[i+1], hist[i]+1)
        linspace_list.extend(temp[0:-1])
    linspace_list.extend(np.linspace(bin_levels[-2], bin_levels[-1], hist[-1]))

    q1, q2, q3 = np.percentile(linspace_list, (25, 50, 75))

    center_prods = [val * cen for val, cen in zip(hist, center_list)]
    center_avg = sum(center_prods)/sum(hist)
    square_diffs = [(val - center_avg) ** 2 for val in center_list]
    sum_square_diffs = sum([n * sqd for n, sqd in zip(hist, square_diffs)])
    var = sum_square_diffs / sum(hist)
    stddev = math.sqrt(var)

    sol_string = "Gjennomsnittet er {:.2f}.\n".format(center_avg)
    sol_string += "Medianen er {:.2f}.\n".format(q2)
    sol_string += "Kvartilbredden er {:.2f}.\n".format(q3-q1)
    sol_string += "Standardavviket er {:.2f}.\n".format(stddev)

    solutions[dictkey] = sol_string

    return return_string
 
 
def four_cummulative_graphs(level="a", seed_string="test", file_name="test.png", sol=False, dictkey=None):
    """Create one bar chart and a figure with four cummulative graphs."""
    bar_chart_file = file_name[:-4] + "-cumbar.png"
    cum_chart_file = file_name[:-4] + "-cumcum.png"
    rnd.seed(seed_string)
     
    bar_list = [round(rnd.gauss(0, 3)) for _ in range(32)]
    values = [bar_list.count(i) for i in range(-10, 10)]
    cum_values = [sum(values[:i+1]) for i in range(len(values))]
    plot_list_bar(values, name=bar_chart_file, start=0, end=20, stepsize=5)
    plot_figure = plt.figure()
 
    ax1 = plot_figure.add_subplot(221)
    ax2 = plot_figure.add_subplot(222)
    ax3 = plot_figure.add_subplot(223)
    ax4 = plot_figure.add_subplot(224)
     
    plots = [
        values[:],
        [value / 2 for value in cum_values],
        [value * 2 for value in cum_values],
        cum_values[:],
        ]
 
    rnd.shuffle(plots)
 
    if sol:
        for num, plot in enumerate(plots):
            if plot == cum_values:
                repl[dictkey] = 'Fasiten er {}.'.format(['a)', 'b)', 'c)', 'd)'][num])
 
    ax1.plot(plots[0])
    ax1.xaxis.set_ticks(np.arange(0, 20, 5))
    ax1.set_title("a)", x=0.1, y=0.8)
    ax2.plot(plots[1])
    ax2.xaxis.set_ticks(np.arange(0, 20, 5))
    ax2.set_title("b)", x=0.1, y=0.8)
    ax3.plot(plots[2])
    ax3.xaxis.set_ticks(np.arange(0, 20, 5))
    ax3.set_title("c)", x=0.1, y=0.8)
    ax4.plot(plots[3])
    ax4.xaxis.set_ticks(np.arange(0, 20, 5))
    ax4.set_title("d)", x=0.1, y=0.8)
 
    plot_figure.savefig("./figurer/" + cum_chart_file)
    plt.close(plot_figure)
    return "[[./figurer/" + bar_chart_file + "]]\n\n[[./figurer/" + cum_chart_file + "]]"
 
 
def pretty_print_list(lst, n=8):
    """Return string of formated values. """
 
    value = ""
    if len(lst) <= n:
        value =  "    ".join(lst)
    else:
        value = ''
        for i in range(0, len(lst), n):
            for j in lst[i:i+n]:
                value += "{0: >8}".format(j)
            value += "\n"
        value = value[:-1]
    return wrap_text(value)
 
 
def wrap_text(string):
    return "#+BEGIN_SRC txt\n" + string +  "\n#+END_SRC"
 
 
def nplusone_average(level="a", seed_string="test", sol=False, dictkey=None):
    """Return string with old average (init), number of students (n) and new average."""
 
    new = rnd.randint(22, 67)
    n = rnd.randint(10, 20)

    init_average = sum(rnd.sample(range(22, 67), n)) / n
    new_average = (init_average * n + new) / (n + 1)

    value = ("På en skole jobber det i utgangspunktet {0} realfagslærere, "
             "som da hadde en gjennomsnittsalder på {1:.2f} år. Det begynner så en "
             "ny lærer på skolen, og brått blir gjennomsnittsalderen {2:.2f} år.\n\n"
             "Hva er aldereden til den nye læreren?").format(n, init_average, new_average)
 
    if sol:
        repl[dictkey] = 'Fasiten er {}.'.format(new)
 
    return value
 
def analyse_points(level="a", length=32,  seed_string="test", sol=False, dictkey=None):
    if level == "a":
        lst = [rnd.randint(1, 20) for _ in range(32)]
 
    elif level == "b":
        lst = [int(rnd.gauss(60, 25)) for _ in range(50)]
        for i, item in enumerate(lst):
            if item < 0:
                lst[i] = 0
            elif item > 100:
                lst[i] = 100
 
 
    if sol and level == "a":
        repl[dictkey] = ("Medianpoengene er {:.2f}.\n"
                                 "Gjennomsnittspoengene er {:.2f}.\n"
                                 "".format(median_list(lst),
                                           average_list(lst)))
    if sol and level == "b":
        repl[dictkey] = analyse_fasit(lst)
 
    return pretty_print_list(lst)
 
def centers(lst):
    """Return a list of length n-1 of centers."""
 
    return [(lst[i] + lst[i+1])/2 for i in range(len(lst) - 1)]
 
def analyse_fasit(lst):
    act_average = average_list(lst)
    bins1 = [0, 20, 40, 60, 80, 95, 100]
    bins2 = [0, 25, 40, 60, 80, 90, 100]
    centers1 = centers(bins1)
    centers2 = centers(bins2)
 
    hist1, bins1_ = np.histogram(lst, bins1)
    hist2, bins2_ = np.histogram(lst, bins2)
 
    return_string = "| Poeng | Antall |\n"
    return_string += "|-|-|\n"
 
    grouped_average = sum([num * points for num, points in zip(hist1, centers1)]) / len(lst)
 
    for i in range(len(bins1)-1):
        return_string += "| [{0}, {1}> | {2} |\n".format(bins1[i], bins1[i+1], hist1[i])
 
    return_string += "\nEksakt gjennomsnitt er {:.2f}.\n".format(act_average)
    return_string += "Gruppert gjennomsnitt er {:.2f}.\n\n".format(grouped_average)
 
    return_string += "Antall 1-ere etter endring er {0}.\n".format(hist2[0])
    return_string += "Antall 6-ere etter endring er {0}.\n".format(hist2[-1])
         
    return return_string
 
 
def analyse_oppgave(level="a", seed_string="test", sol=False, dictkey=None):
    """Returns a string containing the entire text to 'analyseoppgave'."""
    base_a = (
        "I en 2P-klasse ble det gjennomført en prøve, under vises "
        "poengene, som var fra 0 til 20.\n\n"
        "{POINTS}\n\n"
        "- Regn ut gjennomsnitts- og medianpoengene til elevene.\n"
        "- Tegn et diagram som illustrerer dataene"
        "".format(POINTS=analyse_points(level="a", seed_string=seed_string, sol=sol, dictkey=dictkey)))
 
    base_b = (
            "På en skole ble det gjennomført en 2P-prøve for et trinn. Poengene, "
            "som var fra 0 til 100, blei som vist under.\n\n"
            "{POINTS}\n\n"
            "- Grupper dataene i en frekvenstabell hvor gruppene er [0, 20>, [20, 40>, [40, 60>, [60, 80>, [80, 95>, [95, 100].\n"
            "- Finn gjennomsnittspoengene /både/ via det grupperte materialet, og eksakt. Forklar eventuelle forskjeller på de to tallene.\n"
            "- Vis dataene i et histogram.\n\n"
            "Gruppene over viser til en typisk karakterinndeling. En av lærerne på skolen, ønsker å gjøre det litt vanskeligere å bestå,"
            "så alt under 25 poeng gir karakteren 1, men litt lettere å få karakteren 6, nå fra 90 poeng. Hvordan ville det påvirket karakterfordelingen"
            "på trinnet?"
        "".format(POINTS=analyse_points(level="b", seed_string=seed_string, sol=sol, dictkey=dictkey)))
 
    base = (
        "Velg ett av de to alternativene under:\n\n"
        "*Alternativ 1*\n" +
        base_a +
        "*Alternativ 2*\n" +
        base_b
        )

    return base
 
def create_assignment(group="testgroup", student=1, template=template, replacements=None):
 
    name = group + "-" + str(student).zfill(2) + "-" + level
    save_string = template
    image_name = name + ".png"
    page_name = name + ".org"
    page_name_solutions = name + "-fasit.org"

    rnd.seed(name)
 
    save_file = open(page_name, 'w')
    plot_list(name=image_name)

    solution_defaults = {'{SENTRALMAALSOPPGAVEFASIT}' : "",
                         '{FREKVENSTABELLFASITA}' : "",
                         '{FREKVENSTABELLFASITB}' : "",
                         '{NYTTSNITTFASIT}' : "",
                         '{HISTOGRAMFASIT}' : "",
                         '{ANALYSEOPPGAVEFASIT}' : ""}
                       
    solutions = solution_defaults.copy()

    replacements['{FREKVENSTABELLHEADINGB}'] = '| Poeng | Frekvens |'
    replacements['{FREKVENSTABELLHEADINGA}'] = '| Terningkast | Frekvens |'
    replacements['{FREKVENSTABELLDATAA}'] = populate_freq_table_easy(
        seed_string=name, dictkey='{FREKVENSTABELLFASITA }')
    replacements['{FREKVENSTABELLDATAB}'] = populate_freq_table_hard(
        seed_string=name, sol=solutions, dictkey='{FREKVENSTABELLFASITB}')
    replacements['{DEFINISJONSOPPGAVE}'] = draw_words(
        words_easy)
    replacements['{SENTRALMAALSOPPGAVE}'] = pretty_print_list(
        random_list_hard(15, sol=solutions, dictkey='{SENTRALMAALSOPPGAVEFASIT}'))
    replacements['{DIAGRAM}'] = insert_image(
        image_name)
    replacements['{NYTTSNITT}'] = nplusone_average(
        level=level, seed_string=name, sol=solutions, dictkey='{NYTTSNITTFASIT}')
    replacements['{KUMMULATIVMATCH}'] = four_cummulative_graphs(
        level=level, seed_string=name, file_name=image_name, sol=solutions, dictkey='{KUMMULATIVMATCHFASIT}')
    replacements['{ANALYSEOPPGAVE}'] = analyse_oppgave(seed_string=name, sol=solutions, dictkey='{ANALYSEOPPGAVEFASIT}')
    replacements['{HISTOGRAM}'] = plot_histogram(
        histogram_data(level=level, seed_string=name, sol=solutions, dictkey='{HISTOGRAMFASIT}'), file_name=image_name, sol=solutions, dictkey='{HISTOGRAMFASIT}')
    replacements['{HISTOGRAMADD}'] = histogram_addition()
    replacements['{ELEVNUMMER}'] = str(student)
    replacements['{KLASSE}'] = str(group)
    replacements['{LEVEL}'] = str(level)
 
 
    for replacement in replacements:
        save_string = save_string.replace(replacement, replacements[replacement])
     
    save_file.write(save_string)
    save_file.close()
    plt.close('all')
 
    return None
 
print("Generating tasks.")
for group in groups:
    print("Generating for {}".format(group))
    for student in range(1, number_of_students+1):
        print(".", end="")
        repl = {**replacements}
        create_assignment(group, student, replacements=repl)
    print("")
