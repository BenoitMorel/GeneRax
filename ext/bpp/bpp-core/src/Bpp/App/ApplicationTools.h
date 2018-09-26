//
// File: ApplicationTools.h
// Created by: Julien Dutheil
// Created on: Fri Oct 21 16:19 2005
// From old file created on: Sun Dec 14 09:36:26 2003
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

  This software is a computer program whose purpose is to provide basal and 
  utilitary classes. This file belongs to the Bio++ Project.

  This software is governed by the CeCILL  license under French law and
  abiding by the rules of distribution of free software.  You can  use, 
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info". 

  As a counterpart to the access to the source code and  rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty  and the software's author,  the holder of the
  economic rights,  and the successive licensors  have only  limited
  liability. 

  In this respect, the user's attention is drawn to the risks associated
  with loading,  using,  modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean  that it is complicated to manipulate,  and  that  also
  therefore means  that it is reserved for developers  and  experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or 
  data to be ensured and,  more generally, to use and operate it in the 
  same conditions as regards security. 

  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/

#ifndef _APPLICATIONTOOLS_H_
#define _APPLICATIONTOOLS_H_

#include "../Io/FileTools.h"
#include "../Io/OutputStream.h"
#include "../Text/TextTools.h"
#include "../Text/StringTokenizer.h"
#include "../Text/NestedStringTokenizer.h"

#include "../Numeric/Matrix/Matrix.h"

// From the STL:
#include <map>
#include <vector>
#include <iostream>
#include <ctime>

namespace bpp
{

/**
 * @brief This class provides some common tools for developping applications.
 *
 * These functions are designed for helping to parse an option file.
 * 
 * The option files are supposed to follow this simple format:<br>
 * @code
 * parameterName = parameterContent
 * @endcode
 * with one parameter per line.
 *
 * In files, shell comments:
 * @code
 * # my comment line here
 * @endcode
 * C comments:
 * @code
 * / * my comment block here * /
 * @endcode
 * and C++ comments:
 * @code
 * // my comment line here
 * @endcode
 * are allowed, and ignored while parsing.
 *
 * Some methods for displaying information (messages, errors, warnings...) are also provided.
 *
 * Methods dealing with parameters takes as argument a map<string, string> object
 * containing the parameters (names are the keys of the map, and values are... the values of the map!).
 * These map objects may be obtained from the AttributesTools utilitary class.
 */
  class ApplicationTools
  {
  public:
    
    /**
     * @brief The output stream where errors have to be displayed.
     */
    static OutputStream* error;
    /**
     * @brief The output stream where messages have to be displayed.
     */
    static OutputStream* message;
    /**
     * @brief The output stream where warnings have to be displayed.
     */
    static OutputStream* warning;

    /**
     * @brief Timer variable.
     */
    static time_t startTime;

    /**
     * @brief The width of the output terminal (in character).
     */
    static size_t terminalWidth;

    /**
     * @brief The fraction of terminal width dedicated to messages.
     */
    static float terminalSplit;

    /**
     * @brief Tell if the program is interactive (typically run in foreground). Default to yes.
     */
    static bool interactive;
 
    /**
     * @brief Specify the amount of warning to display.
     */
    static int warningLevel;

  public:
    ApplicationTools() {}
    virtual ~ApplicationTools() {}
  
  public:

    /**
     * @brief Tells if a parameter have been specified.
     *
     * @param parameterName The name of the parameter.
     * @param params        The parameter list.
     * @return True is the parameter of specified name is in the list.
     */
    static bool parameterExists(const std::string& parameterName, const std::map<std::string, std::string>& params)
    {
      std::map<std::string, std::string>::const_iterator it=params.find(parameterName);
  
      return (it != params.end() && !TextTools::isEmpty(it->second));
    }


    static bool parameterExists(const std::string& parameterName, std::vector<std::string>& params)
    {
      for (size_t i = 0; i < params.size(); ++i)
        if (params[i] == parameterName)
          return true;

      return false;
    }

    /**
     * @brief Returns a vector of parameter names that match a given pattern.
     *
     * Only "*" wildcard is implemented now.
     *
     * @param pattern The pattern.
     * @param params  The parameter list.
     * @return a vector of matching names.
     */  
    static std::vector<std::string> matchingParameters(const std::string& pattern, std::map<std::string, std::string>& params);

    static std::vector<std::string> matchingParameters(const std::string& pattern, std::vector<std::string>& params);

    /**
     * @brief Get a double parameter.
     *
     * @param parameterName    The name of the corresponding parameter.
     * @param params           The attribute map where options may be found.
     * @param defaultValue     The default value to use if the parameter is not found.
     * @param suffix           A suffix to be applied to the parameter name.
     * @param suffixIsOptional Tell if the suffix is absolutely required.
     * @param warn             Tell if a warning must be sent in case the parameter is not found.
     * @return The corresponding value.
     */
    static double getDoubleParameter(
      const std::string& parameterName,
      std::map<std::string, std::string>& params,
      double defaultValue,
      const std::string& suffix = "",
      bool suffixIsOptional = true,
      int warn = 0);
  
    /**
     * @brief Get an integer parameter.
     *
     * @param parameterName    The name of the corresponding parameter.
     * @param params           The attribute map where options may be found.
     * @param defaultValue     The default value to use if the parameter is not found.
     * @param suffix           A suffix to be applied to the parameter name.
     * @param suffixIsOptional Tell if the suffix is absolutely required.
     * @param warn             Tell if a warning must be sent in case the parameter is not found.
     * @return The corresponding value.
     */
    static int getIntParameter(
      const std::string& parameterName,
      std::map<std::string, std::string>& params,
      int defaultValue,
      const std::string& suffix = "",
      bool suffixIsOptional = true,
      int warn = 0);
  
    /**
     * @brief Get a string parameter.
     *
     * @param parameterName    The name of the corresponding parameter.
     * @param params           The attribute map where options may be found.
     * @param defaultValue     The default value to use if the parameter is not found.
     * @param suffix           A suffix to be applied to the parameter name.
     * @param suffixIsOptional Tell if the suffix is absolutely required.
     * @param warn             Tell if a warning must be sent in case the parameter is not found.
     * @return The corresponding value.
     */
    
    static std::string getStringParameter(
      const std::string& parameterName,
      const std::map<std::string, std::string>& params,
      const std::string& defaultValue,
      const std::string& suffix = "",
      bool suffixIsOptional = true,
      int warn = 0)
    {
      std::string sParam = defaultValue;
      std::map<std::string, std::string>::const_iterator it1=params.find(parameterName + suffix);  
      if (it1 != params.end() && !TextTools::isEmpty(it1->second))
        sParam = it1->second;
      else
      {
        std::map<std::string, std::string>::const_iterator it2=params.find(parameterName);  
        if (suffixIsOptional && it2 != params.end() && !TextTools::isEmpty(it2->second))
          sParam = it2->second;
        else
          if (warn <= warningLevel) {
            displayWarning("Parameter " + parameterName + " not specified. Default used instead: " + defaultValue);
          }
      }
  
      return sParam;
    }
    

    /**
     * @brief Get a boolean parameter.
     *
     * @param parameterName    The name of the corresponding parameter.
     * @param params           The attribute map where options may be found.
     * @param defaultValue     The default value to use if the parameter is not found.
     * @param suffix           A suffix to be applied to the parameter name.
     * @param suffixIsOptional Tell if the suffix is absolutely required.
     * @param warn             Tell if a warning must be sent in case the parameter is not found.
     * @return The corresponding value.
     */
    static bool getBooleanParameter(
      const std::string& parameterName,
      std::map<std::string, std::string>& params,
      bool defaultValue,
      const std::string& suffix = "",
      bool suffixIsOptional = true,
      int warn = 0);

    /**
     * @brief Get a parameter.
     *
     * @param parameterName    The name of the corresponding parameter.
     * @param params           The attribute map where options may be found.
     * @param defaultValue     The default value to use if the parameter is not found.
     * @param suffix           A suffix to be applied to the parameter name.
     * @param suffixIsOptional Tell if the suffix is absolutely required.
     * @param warn             Tell if a warning must be sent in case the parameter is not found.
     * @return The corresponding value.
     */
    template<class T> static T getParameter(
      const std::string& parameterName,
      std::map<std::string, std::string>& params,
      T defaultValue,
      const std::string& suffix = "",
      bool suffixIsOptional = true,
      int warn = 0)
    {
      T tParam = defaultValue;
      if (parameterExists(parameterName + suffix, params))
      {
        tParam = TextTools::to<T>(params[parameterName + suffix]);
      }
      else if (suffixIsOptional && parameterExists(parameterName, params))
      {
        tParam = TextTools::to<T>(params[parameterName]);
      }
      else if (warn <= warningLevel)
      {
        displayWarning("Parameter " + parameterName + suffix + " not specified. Default used instead: " + TextTools::toString(defaultValue));
      }
      return tParam;
    }
  

    /**
     * @brief Get a file path.
     *
     * @param parameter        The name of the corresponding parameter.
     * @param params           The attribute map where options may be found.
     * @param isRequired       Tell if this path is strictly required or is optional
     * (in the first case, if the parameter is not found, the programm will
     * send an error and exit).
     * @param mustExist        Tell if the corresponding file must already exist.
     * @param suffix           A suffix to be applied to the parameter name.
     * @param suffixIsOptional Tell if the suffix is absolutely required.
     * @param defaultPath      Path to use if no argument is provided.
     * @param warn             Tell if a warning must be sent in case the parameter is not found.
     * @throw Exception        If no file path is specified and isRequired is
     *                         true, or the file does not exist and mustExist
     *                         is set to true.
     */
    static std::string getAFilePath(
      const std::string& parameter,
      std::map<std::string, std::string>& params,
      bool isRequired = true,
      bool mustExist = true,
      const std::string& suffix = "",
      bool suffixIsOptional = false,
      const std::string& defaultPath = "none", 
      int warn = 0)
    noexcept(false);

    /**
     * @brief Get a vector.
     *
     * @param parameterName    The name of the corresponding parameter.
     * @param params           The attribute map where options may be found.
     * @param separator        The character used to delimit values.
     * @param defaultValue     The default value to use if the parameter is not found.
     * @param suffix           A suffix to be applied to the parameter name.
     * @param suffixIsOptional Tell if the suffix is absolutely required.
     * @param warn             Tell if a warning must be sent in case the parameter is not found.
     * @return The corresponding value.
     */
    template<class T> static std::vector<T> getVectorParameter(
      const std::string& parameterName,
      std::map<std::string, std::string>& params,
      char separator,
      const std::string& defaultValue,
      const std::string& suffix = "",
      bool suffixIsOptional = true,
      int warn = 0)
    {
      std::string s = getStringParameter(parameterName, params, defaultValue, suffix, suffixIsOptional, warn);
      if (TextTools::isEmpty(s)) return std::vector<T>(0);
      if (s[0] == '(' && s[s.size() - 1] == ')') {
        //This is a delimited vector:
        s = s.substr(1, s.size() - 2);
        if (TextTools::isEmpty(s)) return std::vector<T>(0);
      }
      NestedStringTokenizer st(s, "(", ")", TextTools::toString(separator));
      size_t n = st.numberOfRemainingTokens();
      std::vector<T> v(n);
      for (size_t i = 0; i < n; i++)
      {
        v[i] = TextTools::fromString<T>(st.nextToken());
      }
      return v;
    }

    /**
     * @brief Get a vector.
     *
     * Similar to getVectorParameter, but dedicated to numerical values.
     * It allows the possibility to set range of values, which will be incremented by 1 (like the : operator in R).
     *
     * @param parameterName    The name of the corresponding parameter.
     * @param params           The attribute map where options may be found.
     * @param separator        The character used to delimit values.
     * @param rangeOperator    The character used to delimit ranges (the + 1 operator must be available for T).
     * @param defaultValue     The default value to use if the parameter is not found.
     * @param suffix           A suffix to be applied to the parameter name.
     * @param suffixIsOptional Tell if the suffix is absolutely required.
     * @param warn             Tell if a warning must be sent in case the parameter is not found.
     * @return The corresponding value.
     */
    template<class T> static std::vector<T> getVectorParameter(
      const std::string& parameterName,
      std::map<std::string, std::string>& params,
      char separator,
      char rangeOperator,
      const std::string& defaultValue,
      const std::string& suffix = "",
      bool suffixIsOptional = true,
      bool warn = true)
    {
      std::string s = getStringParameter(parameterName, params, defaultValue, suffix, suffixIsOptional, warn);
      if (s[0] == '(' && s[s.size() - 1] == ')') {
        //This is a delimited vector:
        s = s.substr(1, s.size() - 2);
        if (TextTools::isEmpty(s)) return std::vector<T>(0);
      }
      StringTokenizer st(s, TextTools::toString(separator));
      size_t n = st.numberOfRemainingTokens();
      std::vector<T> v;
      for (size_t i = 0; i < n; i++)
      {
        std::string token = st.nextToken();
        std::string::size_type pos = token.find(rangeOperator);
        if (pos == std::string::npos)
          v.push_back(TextTools::fromString<T>(token));
        else
        {
          T d1 = TextTools::fromString<T>(token.substr(0, pos));
          T d2 = TextTools::fromString<T>(token.substr(pos + 1));
          for (T j = d1; j < d2; j++)
          {
            v.push_back(j);
          }
          v.push_back(d2);
        }
      }
      return v;
    }

    /**
     * @brief Get a RowMatrix. The input is made of embedded
     * parenthesis, such as ((1,2),(3,4)), where the matrix is filled by
     * lines. Here, the matrix would be:
     * \f[
     * \begin{pmatrix}
     * 1 & 2 \\
     * 3 & 4 \\
     * \end{pmatrix}
     * \f]
     *
     * @param parameterName    The name of the corresponding parameter.
     * @param params           The attribute map where options may be found.
     * @param separator        The character used to delimit values.
     * @param defaultValue     The default value to use if the parameter is not found.
     * @param suffix           A suffix to be applied to the parameter name.
     * @param suffixIsOptional Tell if the suffix is absolutely required.
     * @param warn             Tell if a warning must be sent in case the parameter is not found.
     * @return The corresponding value.
     */

    template<class T> static RowMatrix<T> getMatrixParameter(
      const std::string& parameterName,
      std::map<std::string, std::string>& params,
      char separator,
      const std::string& defaultValue,
      const std::string& suffix = "",
      bool suffixIsOptional = true,
      bool warn = true)
    {
      RowMatrix<T> mat;

      std::string s = getStringParameter(parameterName, params, defaultValue, suffix, suffixIsOptional, warn);
      if (TextTools::isEmpty(s)) return RowMatrix<T>(0,0);
      if (s[0] == '(' && s[s.size() - 1] == ')') {
        //This is a delimited vector:
        s = s.substr(1, s.size() - 2);
        if (TextTools::isEmpty(s)) return RowMatrix<T>(0,0);
      }
    
      StringTokenizer st1(s, "()");
    
      while (st1.hasMoreToken())
      {
        std::string si=st1.nextToken();
        StringTokenizer st2(si, TextTools::toString(separator));
        size_t n = st2.numberOfRemainingTokens();

        std::vector<T> v(n);
        for (size_t i = 0; i < n; i++)
        {
          v[i] = TextTools::fromString<T>(st2.nextToken());
        }
      
        if (v.size()!=0)
          mat.addRow(v);
      }
      return mat;
    }


    /**
     * @name Output methods.
     *
     * @{
     */
    
    /**
     * @brief Print a message.
     *
     * @param text The text of the message.
     */
    static void displayMessage(const std::string& text);
    
    /**
     * @brief Print an error message.
     *
     * @param text The text of the message.
     */
    static void displayError(const std::string& text);
    
    /**
     * @brief Print a warning message.
     *
     * @param text The text of the message.
     */
    static void displayWarning(const std::string& text);
    
    /**
     * @brief Print a task message.
     *
     * Display the message and flush the buffer, but do not end the current line.
     *
     * @param text The text of the message.
     * @param eof  Insert a carriage return after displaying the message.
     */
    static void displayTask(const std::string& text, bool eof = false);
    
    /**
     * @brief Print a task ended message.
     *
     * Print "Done." and go to next line.
     */
    static void displayTaskDone();
    
    /**
     * @brief Print a result message.
     *
     * Result will be aligned to 30 character from the begining of the message.
     * ex: text = "Here is what you get:" and result = "THAT" gives
     * "Here is what you get:          THAT".
     *    
     * @param text   The text of the message.
     * @param result The result.
     */
    template<class T>
    static void displayResult(const std::string& text, const T& result)
    {
      displayMessage(TextTools::resizeRight(text, static_cast<size_t>(static_cast<float>(terminalWidth) * terminalSplit - 1), '.') + ": " + TextTools::toString<T>(result));
    }

    /**
     * @brief Print a boolean result message ("yes" or "no").
     *
     * Result will be aligned to 30 character from the begining of the message.
     * @param text   The text of the message.
     * @param result The result.
     */
    static void displayBooleanResult(const std::string& text, bool result)
    {
      displayResult(text, result ? std::string("yes") : std::string("no"));
    }

    /**
     * @brief Display a gauge.
     *
     * Show progress status.
     * @code
     * for(size_t i = 0; i < 1000; i++)
     * {
     *   ApplicationTools::displayGauge(i, 999, '*');
     *   //Perform time consuming task...
     * }
     * @endcode
     * will result in something like:
     * @verbatim
     * [************************************]
     * @endverbatim
     * 
     * @param iter   The current iteration number.
     * @param total  The total number of iteration.
     * @param symbol The character to display in the gauge.
     * @param mes    A message to print before the gauge.
     */
    static void displayGauge(size_t iter, size_t total, char symbol='>', const std::string& mes="");

    /**
     * @brief Display a gauge for unefined amount of iterations.
     *
     * Show progress status.
     * @code
     * for(size_t i = 0; i < 1000; i++)
     * {
     *   ApplicationTools::displayUnlimitedGauge(i);
     *   //Perform time consuming task...
     * }
     * @endcode
     * will result in something like:
     * @verbatim
     * - 1
     * / 2
     * - 3
     * \ 4
     * - 5
     * etc
     * @endverbatim
     * 
     * @param iter   The current iteration number.
     * @param mes    A message to print before the gauge.
     */
    static void displayUnlimitedGauge(size_t iter, const std::string& mes="");


    /** @} */

    /**
     * @brief Starts the timer.
     */
    static void startTimer()
    {
      time(&startTime);
    }

    /**
     * @brief Display the current timer value to the 'message' stream.
     *
     * @param msg Message to display before time.
     */
    static void displayTime(const std::string& msg);
 
    /**
     * @brief Get the current timer value.
     *
     * @return The number of seconds from when timer was started.
     */
    static double getTime();
  };

} //end of namespace bpp.

#endif  //_APPLICATIONTOOLS_H_

