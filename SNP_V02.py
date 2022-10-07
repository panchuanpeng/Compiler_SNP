import os
import tkinter.filedialog
import pandas as pd
import matplotlib.pyplot as plt
from tkinter import ttk
from tkinter import messagebox
from tkinter import *


def full_directory(path, filetype):
	result_dict = {}
	for root1, directory, files in os.walk(path):
		result_list = []
		for file in files:
			if file.endswith(filetype.upper()) or file.endswith(filetype.lower()):
				result_list.append(file)
		result_dict[root1] = result_list

		if len(result_dict) == 0:
			raise Exception('It is an empty directory. No SNP files are existing')
	return result_dict


def get_s2p_skipline(filename):
	count = 1
	with open(filename) as f:
		while count < 20:
			if f.readline()[0] == "#":
				if count == 1:
					if "CMT" in f.readline()[:10]:
						return 5, "\t"
					else:
						raise Exception("it is not a CMT snp file")
					pass
				else:
					if len(f.readline().split(' ')) > 1:
						return count, " "
					elif len(f.readline().split('\t')) > 1:
						return count, "\t"
					else:
						raise Exception("File sep character judgement error")
			else:
				pass
			count += 1
		if count == 20:
			raise Exception("file skipline count error")


def Get_S2P_Column_to_Write(file, column_type, skip_line, sep):
	column = pd.read_csv(file, header=None, sep=sep, skiprows=list(range(skip_line)), usecols=[column_type])
	return column


def S2P_Curve():
	try:
		linetype = combobox.get()
		stype = combobox1.get()
		freq, s11, s21, s12, s22 = 0, 1, 3, 5, 7
		names_dict = full_directory(path1, 's2p')
		for key, VALUES in names_dict.items():
			if len(VALUES) == 0:
				pass
			else:
				plt.figure()
				print(VALUES)
				skip_line, sep = get_s2p_skipline(key + '\\' + VALUES[0])
				for element in VALUES:
					df_freq = Get_S2P_Column_to_Write(key + '\\' + element, freq, skip_line, sep)
					df_db = Get_S2P_Column_to_Write(key + '\\' + element, eval(stype), skip_line, sep)
					df_summary = pd.concat([df_freq, df_db], join='inner', axis=1)
					df_summary.columns = ['freq'] + [element]
					if linetype == "scatter":
						plt.scatter(x=df_summary.get("freq"), y=df_summary.get(element), s=10, label=element)
						plt.legend()
					elif linetype == "line":
						plt.plot(df_summary.get("freq"), df_summary.get(element), label=element)
						plt.legend()
				plt.title(os.path.basename(key) + ": " + stype.upper())
				plt.grid(which='major', axis='both', linestyle='--')
		plt.show()
	except Exception as e:
		messagebox.showerror("Error", message="Error: %s" % e)


def find_file(names_dict, name):
	for key, VALUES in names_dict.items():
		flag, file = 1, ""
		for element in VALUES:
			if element.split(".")[0].lower() == name.lower():
				flag = 0
				file = key + '/' + element
		if flag:
			raise Exception("Not find %s" % name)
		else:
			return file


def S2P_Cal_State():
	try:
		linetype = combobox.get()
		stype = combobox1.get()
		freq, s11, s21, s12, s22 = 0, 1, 3, 5, 7
		names_dict1 = full_directory(path1, 's2p')
		names_dict2 = full_directory(path2, 's2p')
		name_list = ["open", "short", "load", "thru"]
		for name in name_list:
			VALUES = [find_file(names_dict1, name), find_file(names_dict2, name)]
			names = [i.split("/")[-2] + "-" + i.split("/")[-1] for i in VALUES]
			print(VALUES)
			skip_line, sep = get_s2p_skipline(VALUES[0])
			plt.figure()
			for index, element in enumerate(VALUES):
				df_freq = Get_S2P_Column_to_Write(element, freq, skip_line, sep)
				df_db = Get_S2P_Column_to_Write(element, eval(stype), skip_line, sep)
				df_summary = pd.concat([df_freq, df_db], join='inner', axis=1)
				df_summary.columns = ['freq'] + [os.path.split(element)[-1]]
				if linetype == "scatter":
					plt.scatter(x=df_summary.get("freq"), y=df_summary.get(os.path.split(element)[-1]), s=10, label=names[index])
					plt.legend()
				elif linetype == "line":
					plt.plot(df_summary.get("freq"), df_summary.get(os.path.split(element)[-1]), label=names[index])
					plt.legend()
			plt.title(name)
			plt.grid(which='major', axis='both', linestyle='--')
		plt.show()
	except Exception as e:
		messagebox.showerror("Error", message="Error: %s" % e)


def name_for_df(file):
	name_start, name_end = 0, None
	return '_'.join(file.split('\\')[-1].split('_')[name_start:name_end])
	# 需要注意文件命名形式的不一致。根据实际情况调整重命名代码，FAB出来数据为'_'.join(name.split('_')[1:4],LAB数据为'_'.join(name.split('_')[0:4])


def Df_Generation(file, skip_line, sep):
	df_raw = pd.read_csv(file, header=None, sep=sep, skiprows=range(skip_line), skipinitialspace=True)  # 读取该文件的Freq与S21列，组合成df_raw
	df_1 = df_raw.iloc[lambda x: df_raw.index % 3 == 0].reset_index().iloc[:, 1:]
	df_2 = df_raw.iloc[lambda x: df_raw.index % 3 == 1].reset_index().iloc[:, 1:7]
	df_3 = df_raw.iloc[lambda x: df_raw.index % 3 == 2].reset_index().iloc[:, 1:7]
	df_temp_summary = pd.concat([df_1, df_2, df_3], join='inner', axis=1)
	if len(df_1 == len(df_2)):
		return df_temp_summary
	else:
		raise Exception("length of sub dataframe(df_1/2/3) are not aligned")


def get_s3p_skipline(filename):
	count = 1
	with open(filename) as f:
		while count < 20:
			info = f.readline()
			if count == 1:
				test_type = info.split(",")
				if len(test_type) > 1:
					test_type = test_type[1]
			if info[0] == "#":
				if len(info.split(' ')) > 1:
					return test_type, count, " "
				elif len(info.split('\t')) > 1:
					return test_type, count, "\t"
				else:
					raise Exception("File sep character judgement error")
			count += 1
		if count == 20:
			raise Exception("file skipline count error")


def S3P_Curve():
	global num
	try:
		step = 0
		freq, s11, s12, s13, s21, s22, s23, s31, s32, s33 = 0, 1, 3, 5, 7, 9, 11, 13, 15, 17
		linetype = combobox.get()
		stype = combobox2.get()
		names_dict = full_directory(path1, 's3p')
		for key, VALUES in names_dict.items():
			for name in VALUES:
				test_type, skip_line, sep = get_s3p_skipline(key + '\\' + name)
				df_temp = Df_Generation(key + '\\' + name, skip_line, sep)
				if stype == "normal":
					step = 3
					freq, TX_RX, ANT_RX, TX_ANT = freq, s21, s23, s31
					df_TX_RX = df_temp.iloc[:, [freq, TX_RX]]  # TX-RX
					df_ANT_RX = df_temp.iloc[:, [freq, ANT_RX]]  # ANT-RX
					df_TX_ANT = df_temp.iloc[:, [freq, TX_ANT]]  # TX-ANT

					df_TX_RX.columns = ['freq'] + [name_for_df(name)]
					df_ANT_RX.columns = ['freq'] + [name_for_df(name)]
					df_TX_ANT.columns = ['freq'] + [name_for_df(name)]

					plt.figure(num)
					if linetype == "scatter":
						plt.scatter(x=df_TX_ANT.get("freq"), y=df_TX_ANT.get(name_for_df(name)), s=10, label=name_for_df(name))
						plt.legend()
					elif linetype == "line":
						plt.plot(df_TX_ANT.get("freq"), df_TX_ANT.get(name_for_df(name)), label=name_for_df(name))
						plt.legend()
					plt.title(os.path.basename(key) + ': ' + 'TX-ANT')
					plt.grid(which='major', axis='both', linestyle='--')

					plt.figure(num+1)
					if linetype == "scatter":
						plt.scatter(x=df_ANT_RX.get("freq"), y=df_ANT_RX.get(name_for_df(name)), s=10, label=name_for_df(name))
						plt.legend()
					elif linetype == "line":
						plt.plot(df_ANT_RX.get("freq"), df_ANT_RX.get(name_for_df(name)), label=name_for_df(name))
						plt.legend()
					plt.title(os.path.basename(key) + ': ' + 'ANT-RX')
					plt.grid(which='major', axis='both', linestyle='--')

					plt.figure(num + 2)
					if linetype == "scatter":
						plt.scatter(x=df_TX_RX.get("freq"), y=df_TX_RX.get(name_for_df(name)), s=10, label=name_for_df(name))
						plt.legend()
					elif linetype == "line":
						plt.plot(df_TX_RX.get("freq"), df_TX_RX.get(name_for_df(name)), label=name_for_df(name))
						plt.legend()
					plt.title(os.path.basename(key) + ': ' + 'TX-RX')
					plt.grid(which='major', axis='both', linestyle='--')
				else:
					step = 1
					freq, TX_RX = freq, eval(stype)
					df_TX_RX = df_temp.iloc[:, [freq, TX_RX]]  # TX-RX
					df_TX_RX.columns = ['freq'] + [name_for_df(name)]
					plt.figure(num)
					if linetype == "scatter":
						plt.scatter(x=df_TX_RX.get("freq"), y=df_TX_RX.get(name_for_df(name)), s=10, label=name_for_df(name))
						plt.legend()
					elif linetype == "line":
						plt.plot(df_TX_RX.get("freq"), df_TX_RX.get(name_for_df(name)), label=name_for_df(name))
						plt.legend()
					plt.title(os.path.basename(key) + ': ' + stype.upper())
					plt.grid(which='major', axis='both', linestyle='--')
			num = num + step
		plt.show()
	except Exception as e:
		messagebox.showerror("Error", message="Error: %s" % e)


def S3P_Cal_State():
	global num
	try:
		step = 0
		freq, s11, s12, s13, s21, s22, s23, s31, s32, s33 = 0, 1, 3, 5, 7, 9, 11, 13, 15, 17
		linetype = combobox.get()
		stype = combobox2.get()
		names_dict1 = full_directory(path1, 's3p')
		names_dict2 = full_directory(path2, 's3p')
		name_list = ["open", "short", "load", "thru12", "thru13", "thru23"]
		for name in name_list:
			VALUES = [find_file(names_dict1, name), find_file(names_dict2, name)]
			names = [i.split("/")[-2] + "-" + i.split("/")[-1] for i in VALUES]
			print(names)
			for index, element in enumerate(VALUES):
				test_type, skip_line, sep = get_s3p_skipline(element)
				df_temp = Df_Generation(element, skip_line, sep)

				if stype == "normal":
					step = 3
					freq, TX_RX, ANT_RX, TX_ANT = freq, s21, s23, s31
					df_TX_RX = df_temp.iloc[:, [freq, TX_RX]]  # TX-RX
					df_ANT_RX = df_temp.iloc[:, [freq, ANT_RX]]  # ANT-RX
					df_TX_ANT = df_temp.iloc[:, [freq, TX_ANT]]  # TX-ANT

					df_TX_RX.columns = ['freq'] + [names[index]]
					df_ANT_RX.columns = ['freq'] + [names[index]]
					df_TX_ANT.columns = ['freq'] + [names[index]]

					if name == "thru13":
						plt.figure(num)
						if linetype == "scatter":
							plt.scatter(x=df_TX_ANT.get("freq"), y=df_TX_ANT.get(names[index]), s=10, label=names[index])
							plt.legend()
						elif linetype == "line":
							plt.plot(df_TX_ANT.get("freq"), df_TX_ANT.get(names[index]), label=names[index])
							plt.legend()
						plt.title(name + ': ' + 'TX-ANT')
						plt.grid(which='major', axis='both', linestyle='--')
					if name == "thru23":
						plt.figure(num+1)
						if linetype == "scatter":
							plt.scatter(x=df_ANT_RX.get("freq"), y=df_ANT_RX.get(names[index]), s=10,
										label=names[index])
							plt.legend()
						elif linetype == "line":
							plt.plot(df_ANT_RX.get("freq"), df_ANT_RX.get(names[index]), label=names[index])
							plt.legend()
						plt.title(name + ': ' + 'ANT-RX')
						plt.grid(which='major', axis='both', linestyle='--')
					if name == "thru12":
						plt.figure(num+2)
						if linetype == "scatter":
							plt.scatter(x=df_TX_RX.get("freq"), y=df_TX_RX.get(names[index]), s=10, label=names[index])
							plt.legend()
						elif linetype == "line":
							plt.plot(df_TX_RX.get("freq"), df_TX_RX.get(names[index]), label=names[index])
							plt.legend()
						plt.title(name + ': ' + 'TX_RX')
						plt.grid(which='major', axis='both', linestyle='--')
				else:
					step = 1
					freq, TX_RX = freq, eval(stype)
					df_TX_RX = df_temp.iloc[:, [freq, TX_RX]]  # TX-RX
					df_TX_RX.columns = ['freq'] + [names[index]]
					plt.figure(num)
					if linetype == "scatter":
						plt.scatter(x=df_TX_RX.get("freq"), y=df_TX_RX.get(names[index]), s=10, label=names[index])
						plt.legend()
					elif linetype == "line":
						plt.plot(df_TX_RX.get("freq"), df_TX_RX.get(names[index]), label=names[index])
						plt.legend()
					plt.title(name + ': ' + stype.upper())
					plt.grid(which='major', axis='both', linestyle='--')
			num = num + step
		plt.show()
	except Exception as e:
		messagebox.showerror("Error", message="Error: %s" % e)


def get_devive_path1():
	global path1
	path1 = tkinter.filedialog.askdirectory()
	old_path_entry1["text"] = path1


def get_devive_path2():
	global path2
	path2 = tkinter.filedialog.askdirectory()
	old_path_entry2["text"] = path2


if __name__ == "__main__":
	num = 0
	path1, path2 = "", ""
	root = tkinter.Tk()  # 实例化出一个父窗口
	root.title("SNP Editor Designed By PCP")
	screenwidth = root.winfo_screenwidth()
	screenheight = root.winfo_screenheight()
	root.geometry('480x400+%d+%d' % ((screenwidth - 480) / 2, (screenheight - 500) / 2))
	root.resizable(width=False, height=False)
	root.wm_attributes('-topmost', 1)

	old_path_label = Label(root, anchor="w", width=20, font=("", 14), text="Path1: ")
	old_path_label.place(x=60, y=41, anchor='w')
	old_path_entry1 = Label(root, anchor="nw", width=60, font=("", 14), wraplength=300, justify='left')
	old_path_entry1.place(x=120, y=30, anchor='nw')

	old_path_label = Label(root, anchor="w", width=20, font=("", 14), text="Path2: ")
	old_path_label.place(x=60, y=121, anchor='w')
	old_path_entry2 = Label(root, anchor="nw", width=60, font=("", 14), wraplength=300, justify='left')
	old_path_entry2.place(x=120, y=110, anchor='nw')

	value = StringVar()
	values = ['line', 'scatter']
	combobox = ttk.Combobox(master=root, height=10, width=7, state='readonly', cursor='arrow', font=('', 14), textvariable=value, values=values)
	combobox.current(0)
	combobox.place(x=60, y=350, anchor='w')

	value1 = StringVar()
	values = ['s21', 's11', "s12", "s22"]
	combobox1 = ttk.Combobox(master=root, height=10, width=7, state='readonly', cursor='arrow', font=('', 14), textvariable=value1, values=values)
	combobox1.current(0)
	combobox1.place(x=200, y=350, anchor='w')

	value2 = StringVar()
	values = ['normal', 's11', "s12", "s13", 's21', "s22", "s23", 's31', "s32", "s33"]
	combobox2 = ttk.Combobox(master=root, height=10, width=7, state='readonly', cursor='arrow', font=('', 14), textvariable=value2, values=values)
	combobox2.current(0)
	combobox2.place(x=335, y=350, anchor='w')
	# 按钮
	get_path1 = Button(root, text="选择1", bg="lightblue", width=12, command=get_devive_path1)  # 调用内部方法  加()为直接调用
	get_path1.place(x=60, y=220, anchor='w')
	get_path2 = Button(root, text="选择2", bg="lightblue", width=12, command=get_devive_path2)  # 调用内部方法  加()为直接调用
	get_path2.place(x=60, y=270, anchor='w')

	get_button = Button(root, text="S2P_Curve", bg="lightblue", width=12, command=S2P_Curve)  # 调用内部方法  加()为直接调用
	get_button.place(x=200, y=220, anchor='w')
	export_button = Button(root, text="S2P_Calibrate", bg="lightblue", width=12, command=S2P_Cal_State)  # 调用内部方法  加()为直接调用
	export_button.place(x=200, y=270, anchor='w')
	start_button = Button(root, text="S3P_Curve", bg="lightblue", width=12, command=S3P_Curve)  # 调用内部方法  加()为直接调用
	start_button.place(x=335, y=220, anchor='w')
	start_button = Button(root, text="S3P_Calibrate", bg="lightblue", width=12, command=S3P_Cal_State)  # 调用内部方法  加()为直接调用
	start_button.place(x=335, y=270, anchor='w')
	root.mainloop()
