#pragma once

namespace Kratos {
    template <class TDataType>
    class Variable {
        int mKey;
    public:
        Variable(std::string Name) {
            std::hash<std::string> hasher;
            mKey = hasher(Name);
        }

        bool operator==(const Variable& rOther) const {
            return mKey == rOther.mKey;
        }
    };

    static Variable<double> TEMPERATURE("TEMPERATURE");
    static Variable<double> CONDUCTIVITY("CONDUCTIVITY");
    static Variable<double> DENSITY("DENSITY");
}
