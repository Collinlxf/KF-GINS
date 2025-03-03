/*
 * OB_GINS: An Optimization-Based GNSS/INS Integrated Navigation System
 *
 * Copyright (C) 2022 i2Nav Group, Wuhan University
 *
 *     Author : Hailiang Tang
 *    Contact : thl@whu.edu.cn
 *
 * This progrgiengineam is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef FILEBASE_H
#define FILEBASE_H

#include <fstream>
#include <vector>

using std::string;
using std::vector;

class FileBase {
 public:
    /* 静态常量的好处
    1. 共享常量：静态常量在类的所有对象之间共享。这意味着所有对象都可以访问同一个常量，而不是为每个对象单独存储一份常量值。这对于节省内存和保持一致性非常有用。
    2. 类级别访问：静态常量可以在不创建对象的情况下通过类名直接访问。这样可以方便地使用这些常量，而不需要实例化对象。例如，可以通过 FileBase::TEXT 和 FileBase::BINARY 直接访问这些常量。
    3. 编译时常量：静态常量在编译时就确定了值，因此可以用于编译时常量表达式中，比如数组大小等。
     */
    static const int TEXT   = 0;
    static const int BINARY = 1;

    FileBase() = default;
    ~FileBase() {
        if (isOpen())
            filefp_.close();
    }

    void close() {
        filefp_.close();
    }

    bool isOpen() {
        return filefp_.is_open();
    }

    bool isEof() {
        return filefp_.eof();
    }

    std::fstream &fstream() {
        return filefp_;
    }

    int columns() const {
        return columns_;
    }

    /*受保护的成员变量（protected）可以在继承类中直接访问和使用。*/
 protected:
    std::fstream filefp_;
    int filetype_ = TEXT;

    int columns_ = 0;
};

#endif // FILEBASE_H
