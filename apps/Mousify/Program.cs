using Microsoft.Win32;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Mousify
{
    class Program
    {
        static void Main(string[] args)
        {
            var keyName = @"Software\Microsoft\Windows\CurrentVersion\Run";
            using (RegistryKey key = Registry.CurrentUser.OpenSubKey(keyName, true))
            {
                key.SetValue("Mousey", @"c:\bolosky\bin\Mousey.exe");
                key.SetValue("bash", @"c:\windows\system32\bash.exe");
                key.DeleteValue(@"com.squirrel.Teams.Teams");
            }
        }
    }
}
