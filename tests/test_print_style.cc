#include "utils/print_styles.h"

int main() {
  LOG(INFO) << "get " << BOLD("BOLD(\"BOLD\")") << " style string";
  LOG(INFO) << "get " << UNDERLINE("UNDERLINE(\"UNDERLINE\")") << " style string";
  LOG(INFO) << "get " << BLINK("BLINK(\"BLINK\")") << " style string";
  LOG(INFO) << "get " << INVERSE("INVERSE(\"INVERSE\")") << " style string";
  LOG(INFO) << "get " << RED("RED(\"RED\")") << " color string";
  LOG(INFO) << "get " << GREEN("GREEN(\"GREEN\")") << " color string";
  LOG(INFO) << "get " << YELLOW("YELLOW(\"YELLOW\")") << " color string";
  LOG(INFO) << "get " << BLUE("BLUE(\"BLUE\")") << " color string";
  LOG(INFO) << "get " << MAGENTA("MAGENTA(\"MAGENTA\")") << " color string";
  LOG(INFO) << "get " << CYAN("CYAN(\"CYAN\")") << " color string";
  LOG(INFO) << "get " << BRED("BRED(\"BRED\")") << " background color string";
  LOG(INFO) << "get " << BGREEN("BGREEN(\"BGREEN\")") << " background color string";
  LOG(INFO) << "get " << BYELLOW("BYELLOW(\"BYELLOW\")") << " background color string";
  LOG(INFO) << "get " << BBLUE("BBLUE(\"BBLUE\")") << " background color string";
  LOG(INFO) << "get " << BMAGENTA("BMAGENTA(\"BMAGENTA\")") << " background color string";
  LOG(INFO) << "get " << BCYAN("BCYAN(\"BCYAN\")") << " background color string";

  std::string style;
  style = "Bold";
  LogBold("Log" << style) << " style";
  style = "Underline";
  LogUnderline("Log" << style) << " style";
  style = "Blink";
  LogBlink("Log" << style) << " style";
  style = "Inverse";
  LogInverse("Log" << style) << " style";

  style = "Red";
  LogRed("Log" << style) << " color";
  style = "Green";
  LogGreen("Log" << style) << " color";
  style = "Yellow";
  LogYellow("Log" << style) << " color";
  style = "Blue";
  LogBlue("Log" << style) << " color";
  style = "Magenta";
  LogMagenta("Log" << style) << " color";
  style = "Cyan";
  LogCyan("Log" << style) << " color";

  style = "RedBackground";
  LogRedBackground("Log" << style) << " color";
  style = "GreenBackground";
  LogGreenBackground("Log" << style) << " color";
  style = "YellowBackground";
  LogYellowBackground("LogBackground" << style) << " color";
  style = "BlueBackground";
  LogBlueBackground("LogBackground" << style) << " color";
  style = "MagentaBackground";
  LogMagentaBackground("LogBackground" << style) << " color";
  style = "CyanBackground";
  LogCyanBackground("Log" << style) << " color";

  return 0;
}
