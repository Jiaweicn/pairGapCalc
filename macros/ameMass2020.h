//created by DeepSeek online service
//processing AME2020 data

// 核素数据结构
struct Nuclide2020 {
   int A;
   int Z;
   string element;
   double mass_excess;
};

// 解析单行数据
bool parseLine2020(const string& line, Nuclide2020& result) {
   // 跳过控制字符和短行
   if (line.length() < 55) return false;

   // 固定列格式解析
   try {
      // 提取关键字段（基于文件格式描述）
      //int NZ = stoi(line.substr(2, 2));   // 列2-4: NZ
      //int N = stoi(line.substr(6, 3));    // 列5-9: N
      int Z = stoi(line.substr(11, 3));    // 列10-14: Z
      int A = stoi(line.substr(16, 3));   // 列15-19: A
      string el = line.substr(20, 2);     // 列21-23: 元素符号
      string mass_str = line.substr(29, 13); // 列29-42: 质量过剩（keV）

      // 清理元素符号中的空格
      el.erase(remove(el.begin(), el.end(), ' '), el.end());
      if (el.empty()) return false;

      // 处理质量过剩中的特殊字符（如#和*）
      size_t invalid_pos = mass_str.find_first_of("#*");
      if (invalid_pos != string::npos) {
          mass_str = mass_str.substr(0, invalid_pos);
      }
      if (mass_str.empty()) return false;

      // 转换为数值
      double mass_excess = stod(mass_str);

      // 填充结果
      result.A = A;
      result.Z = Z;
      result.element = el;
      result.mass_excess = mass_excess;
      return true;
   } catch (...) {
      return false; // 解析失败
   }
}

// 主函数
bool getNucleusMass2020(short a, short z,double &m){
   const char massFile[60]="/home/cai/prjs/pairGapCalc/ameMass2020.txt";//https://www-nds.iaea.org/amdc/ame2020/mass_1.mas20.txt
   ifstream file(massFile);
   if (!file) {
      cerr << "Error: Cannot open "<<massFile <<"\n";
      return false;
   }

   string line;
   bool in_data_section = false;

   // 跳过文件头，定位到数据区域
   while (getline(file, line)) {
      if (line.find("MASS EXCESS") != string::npos) {
         in_data_section = true;
         getline(file, line); // 跳过列标题行
         break;
      }
   }

   if (!in_data_section) {
      return false;
   }

   // 逐行搜索目标核素
   while (getline(file, line)){
      Nuclide2020 nuclide;
      if (parseLine2020(line, nuclide)) {
         if (nuclide.A == a && nuclide.Z == z) {
             //eleName=nuclide.element 
             double mex=nuclide.mass_excess*1e-3;//MeV
   			 m= a * kamu +mex -z*kMe +(14.4381*pow(z,2.39)+ 1.55468e-6*pow(z,5.35)*1e-6);//formula in AME2020
             return true;
           }
       }
   }
   return false;
}
